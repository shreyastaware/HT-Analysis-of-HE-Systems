tic

Book2 = readtable('Book2.xlsx', 'sheet', 'Input_combined');

PHE_in = Book2.PHE_in;
Water_in = Book2.Water_in;
Water_kgs = Book2.water_kgs;
Prod_gas_kgs = Book2.prod_gas_kgs;
run_hour = Book2.run_hour;

for i = 1: height(Book2) 
%% Properties of Producer gas/ Hot Side at T = 48 deg C

T_prod_gas_in = PHE_in(i);
k_prod_gas = 0.02625;
mu_prod_gas = 1.895e-5;
rho_prod_gas = 1.132;
cp_prod_gas = 1007.4;
nu_prod_gas = mu_prod_gas/rho_prod_gas;

%% Properties of Water/ Cold Side at T = 30 deg C

T_water_in = Water_in(i);
k_water = 0.61550458657174;
mu_water = 7.9734522706505e-4;
rho_water = 995.65204668477;
cp_water = 4180.0202402062;

%% PHE & Plate Dimensions

l_phe_or_plate = 0.340;
w_phe = 0.320;
w_each_plate = 0.005;
h_phe_or_plate = 0.300;
k_aluminium = 205;
no_of_plates = 30;
no_of_gaps = no_of_plates-1;

%% Calculation

m_dot_prod_gas = Prod_gas_kgs(i);
m_dot_water = Water_kgs(i);
C_h = m_dot_prod_gas*cp_prod_gas;
C_c = m_dot_water*cp_water;
        
if C_h > C_c
    C_max = C_h;
    C_min = C_c;
else
    C_max = C_c;
    C_min = C_h;
end

A_crossection_plate_w = h_phe_or_plate * w_each_plate;
P_w = 2*(h_phe_or_plate + w_each_plate);
Hydraulic_Dia_w = (4*A_crossection_plate_w) / P_w;

thk_p = ((w_phe - w_each_plate*no_of_plates)/no_of_gaps);
A_crssection_p = l_phe_or_plate*thk_p;
P_pgas = 2*(l_phe_or_plate + thk_p);
Hydraulic_Dia_p = (4*A_crssection_p) / P_pgas;

C_r = C_min/C_max;

%% Water Side Correlations for finding Nusselt No.

v_mean_w = m_dot_water/(rho_water*A_crossection_plate_w*no_of_plates);
Re_water = (rho_water*v_mean_w*Hydraulic_Dia_w)/mu_water;
Pr_water = (mu_water*cp_water)/k_water;

if Re_water <= 2300
    Nu_w = 4.36;
    
    f_water = 96/ Re_water;
    delta_p_water = (f_water * l_phe_or_plate * rho_water * (v_mean_w^2))/(2 * Hydraulic_Dia_w);

else
    if Re_water <= 5e6 %&& Re_water >= 3000
        f_water = 1/((0.790*log(Re_water) - 1.64)^2);
        %f_water = 1/ ((100 * Re_water)^0.25); % For smooth pipes
        Nu_w = ( (f_water/8) * (Re_water - 1000) * Pr_water )/( 1 + ( 12.7 * ((f_water/8) ^ 0.5) * ( (Pr_water^(2/3)) - 1 ) ) );
    else
        syms f_water;
        func = (1 / sqrt(f_water)) - 2 * log( 2.51 / (Re_water * sqrt(f_water) ) ) == 1;
        f_water = solve(func, f_water);
        %f_water = 1/(0.790*log(Re_water) - 1.64);
        %f_water = 1/ ((100 * Re_water)^0.25); % For smooth pipes
        Nu_w = 0.023*(Re_water^0.8)*(Pr_water^0.4);
    end
    
    delta_p_water = (f_water * l_phe_or_plate * rho_water * (v_mean_w^2))/(2 * Hydraulic_Dia_w);
    
end

h_water = k_water*Nu_w/Hydraulic_Dia_w;

%% Gas Side Correlations for finding Nusselt No.

v_mean_p = m_dot_prod_gas/(no_of_gaps*rho_prod_gas*A_crssection_p); 
Re_g = (rho_prod_gas*v_mean_p*Hydraulic_Dia_p)/mu_prod_gas;
Pr_g = (mu_prod_gas*cp_prod_gas)/k_prod_gas;

if Re_g <= 2300
    Nu_g = 4.36;
    
    f_prod_gas = 96/ Re_g;
    delta_p_prod_gas = (f_prod_gas * h_phe_or_plate * rho_prod_gas * (v_mean_p^2))/(2 * Hydraulic_Dia_p);
else
    if Re_g <= 5e6 %&& Re_g >= 3000
        
        f_prod_gas = 1/((0.790*log(Re_g) - 1.64)^2);
        Nu_g = ( (f_prod_gas/8) * (Re_g - 1000) * Pr_g )/( 1 + ( 12.7 * ((f_prod_gas/8) ^ 0.5) * ( (Pr_g^(2/3)) - 1 ) ) );
    else
        syms f_prod_gas;
        eqn = (1 / sqrt(f_prod_gas)) - 2 * log( 2.51 / (Re_g * sqrt(f_prod_gas) ) ) == 1;
        f_prod_gas = solve(eqn, f_prod_gas);
        %f_prod_gas = 1/(0.790*log(Re_g) - 1.64);
        %f_prod_gas = 1/ ((100 * Re_g)^0.25); % For smooth pipes 
        Nu_g = 0.023*(Re_g^0.8)*(Pr_g^0.3);
    end
    
    delta_p_prod_gas = (f_prod_gas * h_phe_or_plate * rho_prod_gas * (v_mean_p^2))/(2 * Hydraulic_Dia_p);
    
end

h_prod_gas = k_prod_gas*Nu_g/Hydraulic_Dia_p;

%% Overall HT Coeff & Eff. NTU Corr

U_t = 1/( (1/h_water) + ((2*0.0001*l_phe_or_plate*h_phe_or_plate)/(k_aluminium*l_phe_or_plate*h_phe_or_plate)) + ((2*0.0001*l_phe_or_plate*h_phe_or_plate)/(k_aluminium*w_each_plate*l_phe_or_plate)) + (1/h_prod_gas) );
A_t = 28*2*((l_phe_or_plate * h_phe_or_plate) + (l_phe_or_plate*w_each_plate)) + 2*(2*w_each_plate*l_phe_or_plate + (l_phe_or_plate * h_phe_or_plate));

NTU = U_t*A_t/C_min;

%Eff = 1 - exp( (-1/C_r) * ( 1 - exp(- C_r * NTU ) ) ) % C_min = mixed and C_max = unmixed
Eff = (1/C_r) * (1 - exp(-C_r * (1 - exp(-NTU)) ) ); % C_min = unmixed and C_max = mixed

q_max = C_min * (T_prod_gas_in - T_water_in);

T_prod_gas_out = T_prod_gas_in - Eff * (T_prod_gas_in - T_water_in);
T_water_out = T_water_in + C_r * (T_prod_gas_in - T_prod_gas_out);

q = Eff * C_min * (T_prod_gas_in - T_water_in);

mdot_gas_1 = q/(cp_prod_gas*(T_prod_gas_in - T_prod_gas_out));

%% Liquid Condensation and RH Calculation

Pgst1 = XSteam('psat_T', T_prod_gas_in);
phi1 = 1;
Pv1 = phi1*Pgst1;
PT2 = (70*10/136 + 710)/760;
PT1 = PT2 + delta_p_prod_gas*1e-5;

Pa1 = PT1 - Pv1;
w1 = 0.88*Pv1/Pa1;
m_dot_g1 = m_dot_prod_gas/(1+w1); 
m_dot_v1 = w1*m_dot_g1;
m_dot_g2 = m_dot_g1;

h_vap_T1 = XSteam('hV_T', T_prod_gas_in);
h_vap_T2 = XSteam('hV_T', T_prod_gas_out);
h_l2 = XSteam('hL_T', T_prod_gas_out);

m_dot_liq = (q - m_dot_g1*cp_prod_gas*(T_prod_gas_in - T_prod_gas_out) - m_dot_v1*1000*(h_vap_T1 - h_vap_T2))/ ((h_vap_T2 - h_l2)*1000);
m_dot_v2 = m_dot_v1 - m_dot_liq;
w2 = m_dot_v2 / m_dot_g2;

Pv2 = w2*PT2/(0.88 + w2);
Pgst2  = XSteam('psat_T', T_prod_gas_out);
%phi2 = Pv2/Pgst2;

if Pgst2 < Pv2
    phi2 = 1;
    Pv2 = phi2*Pgst2;
    Pg2 = PT2 -Pv2;
    w2 = 0.88*Pv2/Pg2;
    m_dot_v2 = w2*m_dot_g2;
    m_dot_liq = m_dot_v1 - m_dot_v2;
    
    %while false
       %if PHE_in >= 50
       %     drainage = m_dot_liq*Eff*run_hour(i)*3600*1000/rho_water;
       %else
       %     drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
       % end 
    %end
    %drainage = m_dot_liq*Eff*run_hour(i)*3600*1000/rho_water;
    drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    
else
    xvap2 =  Pv2 / PT2;
    xsat2 = Pgst2 / PT2;
    phi2 = xvap2/xsat2;
    Pv2 = phi2*Pgst2;
    Pg2 = PT2 -Pv2;
    w2 = 0.88*Pv2/Pg2;
    m_dot_v2 = w2*m_dot_g2;
    m_dot_liq = m_dot_v1 - m_dot_v2;
    
    %while false
    %    if PHE_in >= 50
    %        drainage = m_dot_liq*Eff*run_hour(i)*3600*1000/rho_water;
    %    else
    %        drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    %    end
    %end
    %drainage = m_dot_liq*Eff*run_hour(i)*3600*1000/rho_water;
    drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    
end

Book2.massflwnew(i) = mdot_gas_1/m_dot_prod_gas;
Book2.T_prod_gas_out(i) = T_prod_gas_out;
Book2.T_water_out(i) = T_water_out;
Book2.q(i) = q;
Book2.q_max(i) = q_max;
Book2.Eff(i) = Eff;
Book2.Delta_p_water(i) = delta_p_water;
Book2.Delta_p_prod_gas(i) = delta_p_prod_gas;
Book2.Re_water(i) = Re_water;
Book2.Pr_water(i) = Pr_water;
Book2.Re_prod_gas(i) = Re_g;
Book2.Pr_prod_gas(i) = Pr_g;
Book2.m_dot_condensed(i) = m_dot_liq;
Book2.RH_fraction(i) = phi2;
Book2.DinL(i) = drainage;

end

%Book2(:,:)
toc
