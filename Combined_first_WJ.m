tic

Book2 = readtable('Book2.xlsx', 'sheet', 'Input_combined');

WJ2_in = Book2.PHE_in;
Water_in = Book2.Water_in;
Water_kgs = Book2.water_kgs;
Prod_gas_kgs = Book2.prod_gas_kgs;
run_hour = Book2.run_hour;

for i = 1: height(Book2)
%% Properties of Producer gas/ Hot Side at T = 48 deg C

T_prod_gas_in = WJ2_in(i);
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

%% Multi-pipe Counterflow/ Parallel Flow Heat Exchanger

l_each_tube = 1.116;
D1 = 0.200;
do = 0.100;
thk = 0.005;
di = do - thk;
k_aluminium = 205;
no_of_pipes = 1;

%% Calculation

A_h = no_of_pipes * pi * di * l_each_tube;
A_c = no_of_pipes * pi * do * l_each_tube;
A =  (A_h + A_c)/2;

A_cs_pipe = (pi*(di^2))/4;
A_gas_side = ((no_of_pipes*pi*(do^2))/4);
A_cs_annular_reg = ((pi*(D1^2))/4)- A_gas_side;
P_annular_reg = pi*D1 + no_of_pipes*pi*do;
Hydraulic_Dia_w = 4*A_cs_annular_reg/P_annular_reg;

%m_dot_prod_gas1 = Prod_gas_kgs(i);
%m_dot_prod_gas = 0.5*(m_dot_prod_gas1 + (m_dot_prod_gas1^2 + 2*9.81*l_each_tube*(rho_prod_gas^2)*(A_gas_side^2))^0.5);
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

C_r = C_min/C_max;

%% Water Side Correlations for finding Nusselt No.

v_mean_w = m_dot_water/(rho_water*A_cs_annular_reg);
Re_water = (rho_water*v_mean_w*Hydraulic_Dia_w)/mu_water;
Pr_water = (mu_water*cp_water)/k_water;

if Re_water <= 2300
    Nu_w = 4.36;
    
    f_water = 64/ Re_water;
    delta_p_water = (f_water * l_each_tube * rho_water * (v_mean_w^2))/(2 * Hydraulic_Dia_w) + (rho_water/2) * (v_mean_w^2)*(1 + (1 - pi*0.036*0.036/(4*0.170*1.08  ))^2);
else
    if Re_water <= 5e6 %&& Re_water >= 3000
        f_water = 0.316/ (Re_water^0.25); % For smooth pipes
        f_water = 1/((0.790*log(Re_water) - 1.64)^2);
        Nu_w = ( (f_water/8) * (Re_water - 1000) * Pr_water )/( 1 + ( 12.7 * ((f_water/8) ^ 0.5) * ( (Pr_water^(2/3)) - 1 ) ) );
    else
        syms f_water;
        func = (1 / sqrt(f_water)) - 2 * log10( 2.51 / (Re_water * sqrt(f_water) ) ) == 1;
        f_water = solve(func, f_water);
        %f_water = 1/ ((100 * Re_water)^0.25); % For smooth pipes
        Nu_w = 0.023*(Re_water^0.8)*(Pr_water^0.4);
    end

    delta_p_water = (f_water * l_each_tube * rho_water * (v_mean_w^2))/(2 * Hydraulic_Dia_w) + (rho_water/2) * (v_mean_w^2)*(1 + (1 - pi*0.036*0.036/(4*0.170*1.08  ))^2);
    %delta_p_water = (((f_water * l_each_tube * rho_water * (v_mean_w^2))/(2 * Hydraulic_Dia_w)) - (rho_water * 9.81 * l_each_tube));
end

h_water = k_water*Nu_w/Hydraulic_Dia_w;

%% Gas Side Correlations for finding Nusselt No.

%v_mean_bef = m_dot_prod_gas/(rho_prod_gas*pi*0.170*0.170/4);
v_mean_p = m_dot_prod_gas/(no_of_pipes*rho_prod_gas*A_cs_pipe); 
Re_g = (rho_prod_gas*v_mean_p*di)/mu_prod_gas;
Pr_g = (mu_prod_gas*cp_prod_gas)/k_prod_gas;

if Re_g <= 2300
    Nu_g = 4.36;
    
    f_prod_gas = 64/ Re_g;
    delta_p_prod_gas = ((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2* di)) + (1.5*rho_prod_gas * (v_mean_p^2))/2;
    %delta_p_prod_gas = (((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2 * di)) - (rho_prod_gas * 9.81 * l_each_tube));
else
    if Re_g <= 5e6 %&& Re_g >= 3000
        f_prod_gas = 1/((0.790*log(Re_g) - 1.64)^2);
        f_prod_gas = 1/ ((100 * Re_g)^0.25); % For smooth pipes
        Nu_g = ( (f_prod_gas/8) * (Re_g - 1000) * Pr_g )/( 1 + ( 12.7 * ((f_prod_gas/8) ^ 0.5) * ( (Pr_g^(2/3)) - 1 ) ) );
    else
        syms f_prod_gas;
        eqn = (1 / sqrt(f_prod_gas)) - 2 * log10( 2.51 / (Re_g * sqrt(f_prod_gas) ) ) == 1;
        f_prod_gas = solve(eqn, f_prod_gas);
        %f_prod_gas = 1/ ((100 * Re_g)^0.25); % For smooth pipes
        Nu_g = 0.023*(Re_g^0.8)*(Pr_g^0.3);
    end
    
    delta_p_prod_gas = ((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2* di)) + (1.5*rho_prod_gas * (v_mean_p^2))/2;
    %delta_p_prod_gas = (((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2 * di)) - (rho_prod_gas * 9.81 * l_each_tube));
end

h_prod_gas = k_prod_gas*Nu_g/di;

%% Overall HT Coeff & Eff. NTU Corr

U_t = 1 /( A * ((1/(h_water*A_c)) + ((log(do/di))/(2*no_of_pipes*pi*k_aluminium*l_each_tube)) + (1/(h_prod_gas*A_h)) ) );
NTU = U_t*A/C_min;

Eff_PF = (1 - exp(-NTU * (1 + C_r)) )/(1 + C_r); % Parallel Flow HE

q_max = C_min * (T_prod_gas_in - T_water_in);

T_prod_gas_out_PF = T_prod_gas_in - Eff_PF * (T_prod_gas_in - T_water_in);
T_water_out_PF = T_water_in + C_r * (T_prod_gas_in - T_prod_gas_out_PF);

q_PF = Eff_PF * C_min * (T_prod_gas_in - T_water_in);

Book2.T_prod_gas_out_PF(i) = T_prod_gas_out_PF;
Book2.T_water_out_PF(i) = T_water_out_PF;
Book2.q_PF(i) = q_PF;
Book2.q_max_PF(i) = q_max;
Book2.Eff_PF(i) = Eff_PF;
Book2.Delta_p_water_PF(i) = delta_p_water;
Book2.Delta_p_prod_gas_PF(i) = delta_p_prod_gas;

%% Liquid Condensation and RH Calculation

Pgst1 = XSteam('psat_T', T_prod_gas_in);
phi1 = 1;
Pv1 = phi1*Pgst1;
PT1 = (70*10/136 + 710)/760 + 2*delta_p_prod_gas*1e-5;
delta_p_prod_gas = delta_p_prod_gas + (rho_prod_gas * 9.81 * l_each_tube);
PT2 = (70*10/136 + 710)/760 + delta_p_prod_gas*1e-5;

Pg1 = PT1 - Pv1;
w1 = 0.88*Pv1/Pg1;
m_dot_g1 = m_dot_prod_gas/(1+w1);
m_dot_v1 = w1*m_dot_g1;
m_dot_g2 = m_dot_g1;

h_vap_T1 = XSteam('hV_T', T_prod_gas_in);
h_vap_T2 = XSteam('hV_T', T_prod_gas_out_PF);
h_l2 = XSteam('hL_T', T_prod_gas_out_PF);

m_dot_liq = (q_PF - m_dot_g1*cp_prod_gas*(T_prod_gas_in - T_prod_gas_out_PF) - m_dot_v1*1000*(h_vap_T1 - h_vap_T2))/ ((h_vap_T2 - h_l2)*1000);
m_dot_v2 = m_dot_v1 - m_dot_liq;
w2 = m_dot_v2 / m_dot_g2;

Pv2 = w2*PT2/(0.88 + w2);
Pgst2 = XSteam('psat_T', T_prod_gas_out_PF);

if Pgst2 < Pv2
    phi2 = 1;
    Pv2 = phi2*Pgst2;
    Pg2 = PT2 -Pv2;
    w2 = 0.88*Pv2/Pg2;
    m_dot_v2 = w2*m_dot_g2;
    m_dot_liq = m_dot_v1 - m_dot_v2;
    m_dot_prod_gas_out = m_dot_v2 + m_dot_g2;
    
    %if WJ2_in(i) >= 50
    %    drainage = m_dot_liq*Eff_PF*run_hour(i)*3600*1000/rho_water;
    %else
    %    drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    %end
    %drainage = m_dot_liq*Eff_PF*run_hour(i)*3600*1000/rho_water;
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
    m_dot_prod_gas_out = m_dot_v2 + m_dot_g2;
    
    %if WJ2_in(i) >= 50
    %    drainage = m_dot_liq*Eff_PF*run_hour(i)*3600*1000/rho_water;
    %else
    %    drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    %end
    %drainage = m_dot_liq*Eff_PF*run_hour(i)*3600*1000/rho_water;
    drainage = m_dot_liq*run_hour(i)*3600*1000/rho_water;
end

Book2.m_dot_condensed(i) = m_dot_liq;
Book2.RH_fraction(i) = phi2;
%Book2.Vap_press_inlet(i) = Pv1;
%Book2.Vap_press_outlet(i) = Pv2;
%Book2.Delta_p_prod_gas(i) = delta_p_prod_gas;
Book2.DinL(i) = drainage;
Book2.m_dot_prod_gas_out(i) = m_dot_prod_gas_out;
Book2.PT1(i) = PT1;
Book2.PT2(i) = PT2;

%% Properties of Producer gas/ Hot Side at T = 48 deg C

T_prod_gas_in = T_prod_gas_out_PF;
T_water_in = Water_in(i);

%% Calculation

%m_dot_prod_gas2 = m_dot_prod_gas_out;
%m_dot_prod_gas = 0.5*(m_dot_prod_gas2 + (m_dot_prod_gas2^2 - 2*9.81*l_each_tube*(rho_prod_gas^2)*(A_gas_side^2))^0.5);
m_dot_prod_gas = m_dot_prod_gas_out;
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

C_r = C_min/C_max;

%% Gas Side Correlations for finding Nusselt No.

%v_mean_bef = m_dot_prod_gas/(rho_prod_gas*pi*0.170*0.170/4);
v_mean_p = m_dot_prod_gas/(no_of_pipes*rho_prod_gas*A_cs_pipe);
Re_g = (rho_prod_gas*v_mean_p*di)/mu_prod_gas;
Pr_g = (mu_prod_gas*cp_prod_gas)/k_prod_gas;

if Re_g <= 2300
    Nu_g = 4.36;
    
    f_prod_gas = 64/ Re_g;
    delta_p_prod_gas = ((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2* di)) + (1.5*rho_prod_gas * (v_mean_p^2))/2;
    %delta_p_prod_gas = (((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2 * di)) - (rho_prod_gas * 9.81 * l_each_tube));
else
    if Re_g <= 5e6 %&& Re_g >= 3000
        f_prod_gas = 1/((0.790*log(Re_g) - 1.64)^2);
        f_prod_gas = 1/ ((100 * Re_g)^0.25); % For smooth pipes
        Nu_g = ( (f_prod_gas/8) * (Re_g - 1000) * Pr_g )/( 1 + ( 12.7 * ((f_prod_gas/8) ^ 0.5) * ( (Pr_g^(2/3)) - 1 ) ) );
    else
        syms f_prod_gas;
        eqn = (1 / sqrt(f_prod_gas)) - 2 * log10( 2.51 / (Re_g * sqrt(f_prod_gas) ) ) == 1;
        f_prod_gas = solve(eqn, f_prod_gas);
        %f_prod_gas = 1/ ((100 * Re_g)^0.25); % For smooth pipes
        Nu_g = 0.023*(Re_g^0.8)*(Pr_g^0.3);
    end
    
    delta_p_prod_gas = ((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2* di)) + (1.5*rho_prod_gas * (v_mean_p^2))/2;
    %delta_p_prod_gas = (((f_prod_gas * l_each_tube * rho_prod_gas * (v_mean_p^2))/(2 * di)) - (rho_prod_gas * 9.81 * l_each_tube));
end

h_prod_gas = k_prod_gas*Nu_g/di;

%% Overall HT Coeff & Eff. NTU Corr

U_t = 1 /( A * ((1/(h_water*A_c)) + ((log(do/di))/(2*no_of_pipes*pi*k_aluminium*l_each_tube)) + (1/(h_prod_gas*A_h)) ) );
NTU = U_t*A/C_min;

q_max = C_min * (T_prod_gas_in - T_water_in);

Book2.q_max_CF(i) = q_max;
Book2.Delta_p_water_CF(i) = delta_p_water;
Book2.Delta_p_prod_gas_CF(i) = delta_p_prod_gas;

Eff_CF = (1 - exp(-NTU * (1 - C_r)) )/(1 - (C_r * exp( -NTU * (1 - C_r) ) )); % Counter Flow HE

T_prod_gas_out_CF = T_prod_gas_in - Eff_CF * (T_prod_gas_in - T_water_in);
T_water_out_CF = T_water_in + C_r * (T_prod_gas_in - T_prod_gas_out_CF);

q_CF = Eff_CF * C_min * (T_prod_gas_in - T_water_in);

Book2.T_prod_gas_out_CF(i) = T_prod_gas_out_CF;
Book2.T_water_out_CF(i) = T_water_out_CF;
Book2.q_CF(i) = q_CF;
Book2.Eff_CF(i) = Eff_CF;

%% Liquid Condensation and RH Calculation

Pgst1 = XSteam('psat_T', T_prod_gas_in);
phi1 = phi2;
Pv1 = phi1*Pgst1;
delta_p_prod_gas = delta_p_prod_gas + (rho_prod_gas * 9.81 * l_each_tube);
PT2 = (70*10/136 + 710)/760;
PT1 = PT2 + delta_p_prod_gas*1e-5;

Pg1 = PT1 - Pv1;
w1 = 0.88*Pv1/Pg1;
m_dot_g1 = m_dot_prod_gas/(1+w1);
m_dot_v1 = w1*m_dot_g1;
m_dot_g2 = m_dot_g1;

h_vap_T1 = XSteam('hV_T', T_prod_gas_in);
h_vap_T2 = XSteam('hV_T', T_prod_gas_out_CF);
h_l2 = XSteam('hL_T', T_prod_gas_out_CF);

m_dot_liq = (q_CF - m_dot_g1*cp_prod_gas*(T_prod_gas_in - T_prod_gas_out_CF) - m_dot_v1*1000*(h_vap_T1 - h_vap_T2))/ ((h_vap_T2 - h_l2)*1000);
m_dot_v2 = m_dot_v1 - m_dot_liq;
w2 = m_dot_v2 / m_dot_g2;

Pv2 = w2*PT2/(0.88 + w2);
Pgst2 = XSteam('psat_T', T_prod_gas_out_CF);

if Pgst2 < Pv2
    phi2 = 1;
    Pv2 = phi2*Pgst2;
    Pg2 = PT2 -Pv2;
    w2 = 0.88*Pv2/Pg2;
    m_dot_v2 = w2*m_dot_g2;
    m_dot_liq_1 = m_dot_v1 - m_dot_v2;
    m_dot_prod_gas_out = m_dot_v2 + m_dot_g2;
    
    %if WJ2_in(i) >= 50
    %    drainage_1 = m_dot_liq*Eff_CF*run_hour(i)*3600*1000/rho_water;
    %else
    %    drainage_1 = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    %end
    %drainage_1 = m_dot_liq*Eff_CF*run_hour(i)*3600*1000/rho_water;
    drainage_1 = m_dot_liq*run_hour(i)*3600*1000/rho_water;
else
    xvap2 =  Pv2 / PT2;
    xsat2 = Pgst2 / PT2;
    phi2 = xvap2/xsat2;
    Pv2 = phi2*Pgst2;
    Pg2 = PT2 -Pv2;
    w2 = 0.88*Pv2/Pg2;
    m_dot_v2 = w2*m_dot_g2;
    m_dot_liq_1 = m_dot_v1 - m_dot_v2;
    m_dot_prod_gas_out = m_dot_v2 + m_dot_g2;
    
    %if WJ2_in(i) >= 50
    %    drainage_1 = m_dot_liq*Eff_CF*run_hour(i)*3600*1000/rho_water;
    %else
    %    drainage_1 = m_dot_liq*run_hour(i)*3600*1000/rho_water;
    %end
    %drainage_1 = m_dot_liq*Eff_CF*run_hour(i)*3600*1000/rho_water;
    drainage_1 = m_dot_liq*run_hour(i)*3600*1000/rho_water;
end

dra = drainage_1 + drainage;
Book2.m_dot_condensed_CF(i) = m_dot_liq_1;
Book2.RH_fraction_CF(i) = phi2;
Book2.Delta_p_prod_gas_CF(i) = delta_p_prod_gas;
Book2.DinL_total(i) = dra;
Book2.PT2_Clone(i) = PT1;
Book2.PT3(i) = PT2;

end

%Book2(:,:)
toc