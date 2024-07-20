/////////////////////////////////////////////////////////////
// The first part is just the standard full Hertzian--Mindlin
/////////////////////////////////////////////////////////////

// No need to do any contact force calculation if no contact. And it can happen,
// since we added extra contact margin for adding electrostatic forces before
// physical contacts emerge.
/////////////////////////////////////////////////////////////
// The first part is just the standard full Hertzian--Mindlin
/////////////////////////////////////////////////////////////

// No need to do any contact force calculation if no contact. And it can happen,
// since we added extra contact margin for adding electrostatic forces before
// physical contacts emerge.
float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
const float ks = 2e5;

if (overlapDepth > 0) {
    // Material properties
    // float E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
    {
        // E and nu are associated with each material, so obtain them this way
        float E_A = E[bodyAMatType];
        float nu_A = nu[bodyAMatType];
        float E_B = E[bodyBMatType];
        float nu_B = nu[bodyBMatType];
        matProxy2ContactParam<float>(E_cnt, G_cnt, E_A, nu_A, E_B, nu_B);
        // CoR, mu and Crr are pair-wise, so obtain them this way
        CoR_cnt = CoR[bodyAMatType][bodyBMatType];
        mu_cnt = mu[bodyAMatType][bodyBMatType];
        Crr_cnt = Crr[bodyAMatType][bodyBMatType];
    }

    float3 rotVelCPA, rotVelCPB;
        // We also need the relative velocity between A and B in global frame to use in the damping terms
        // To get that, we need contact points' rotational velocity in GLOBAL frame
        // This is local rotational velocity (the portion of linear vel contributed by rotation)
        rotVelCPA = cross(ARotVel, locCPA);
        rotVelCPB = cross(BRotVel, locCPB);
        // This is mapping from local rotational velocity to global
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPA.x, rotVelCPA.y, rotVelCPA.z, AOriQ.w, AOriQ.x, AOriQ.y,
                                                AOriQ.z);
        applyOriQToVector3<float, deme::oriQ_t>(rotVelCPB.x, rotVelCPB.y, rotVelCPB.z, BOriQ.w, BOriQ.x, BOriQ.y,
                                                BOriQ.z);

    // A few re-usables
    float mass_eff, sqrt_Rd, beta;
    float3 vrel_tan;
    float3 delta_tan = make_float3(delta_tan_x, delta_tan_y, delta_tan_z);

    // Normal force part
        // The (total) relative linear velocity of A relative to B
        const float3 velB2A = (ALinVel + rotVelCPA) - (BLinVel + rotVelCPB);
        const float projection = dot(velB2A, B2A);
        vrel_tan = velB2A - projection * B2A;

        // Now we already have sufficient info to update contact history
        {
            delta_tan += ts * vrel_tan;
            const float disp_proj = dot(delta_tan, B2A);
            delta_tan -= disp_proj * B2A;
            delta_time += ts;
        }

        mass_eff = (AOwnerMass * BOwnerMass) / (AOwnerMass + BOwnerMass);
        sqrt_Rd = sqrt(overlapDepth * (ARadius * BRadius) / (ARadius + BRadius));
        const float Sn = 2. * E_cnt * sqrt_Rd;

        const float loge = (CoR_cnt < DEME_TINY_FLOAT) ? log(DEME_TINY_FLOAT) : log(CoR_cnt);
        beta = loge / sqrt(loge * loge + deme::PI_SQUARED);

        const float k_n = deme::TWO_OVER_THREE * Sn;
        const float gamma_n = deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta * sqrt(Sn * mass_eff);

        force += (k_n * overlapDepth + gamma_n * projection) * B2A;

    // Rolling resistance part

    if (mu_cnt > 0.0) {
        const float kt = 8. * G_cnt * sqrt_Rd;
        const float gt = -deme::TWO_TIMES_SQRT_FIVE_OVER_SIX * beta * sqrt(mass_eff * kt);
        float3 tangent_force = -kt * delta_tan - gt * vrel_tan;
        const float ft = length(tangent_force);
        if (ft > DEME_TINY_FLOAT) {
            // Reverse-engineer to get tangential displacement
            const float ft_max = length(force) * mu_cnt;
            if (ft > ft_max) {
                tangent_force = (ft_max / ft) * tangent_force;
                delta_tan = (tangent_force + gt * vrel_tan) / (-kt);
            }
        } else {
            tangent_force = make_float3(0, 0, 0);
        }
        // Use force to collect tangent_force
        force += tangent_force;
    }

    // Finally, make sure we update those wildcards (in this case, contact history)
    delta_tan_x = delta_tan.x;
    delta_tan_y = delta_tan.y;
    delta_tan_z = delta_tan.z;



    // force magnitude
    double force_mag =  k_n * overlapDepth + gamma_n * projection;

    if (force_mag < 0) {
        force_mag = 0;
    }

    int curr_step = (int)(time / ts);

    if (curr_step % 2000 == 0) {
        // radius contact
        double radius_eff = (ARadius * BRadius) / (ARadius + BRadius);
        double radius_contact = powf( 2.f * radius_eff / E_cnt * force_mag, 1.0/3.0);

        // temperature
        double T_j = Temp_B[BGeo];
        double T_i = Temp_A[AGeo];
        double Q_ij = 2. * ks * radius_contact *  (T_j - T_i);
        atomicAdd(Q_A + AGeo, Q_ij);
        atomicAdd(Q_B + BGeo, -Q_ij);

    }




} else {
    // This is to be more rigorous. If in fact no physical contact, then contact wildcards (such as contact history)
    // should be cleared (they will also be automatically cleared if the contact is no longer detected).
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
    
    }

    // this is where particle fluid particle heat exchange model should be added
    int curr_step = (int) (time / ts);
    if (curr_step % 2000 == 0) {
        double distances[11] = {-0.200, -0.160, -0.120, -0.080, -0.040, 0.000, 0.040, 0.080, 0.120, 0.16, 0.200};
        double volumes[11] = {0.150, 0.167, 0.182, 0.197, 0.217, 0.409, 0.077, 0.046, 0.027, 0.012, 0.000};

        double interval = 0.04;
        double start_distance = -0.2;

        // calculate index
        double ratio = overlapDepth / (ARadius + BRadius);
        int i = (int) ((ratio - start_distance) / interval);

        // Check if the distance is within bounds
        if (i > 0 && i < 11 - 1) {

        double t = (ratio - distances[i]) / (distances[i + 1] - distances[i]);
        double volume = volumes[i] * (1 - t) + volumes[i + 1] * t;


        // temperature
        double T_j = Temp_B[BGeo];
        double T_i = Temp_A[AGeo];
        double Q_ij = ks * volume *  (T_j - T_i);
        atomicAdd(Q_A + AGeo, Q_ij);
        atomicAdd(Q_B + BGeo, -Q_ij);
        }
    }
