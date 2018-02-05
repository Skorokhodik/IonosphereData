function [Vz_interpolated]=Naiquist_theorem(Vz, L, Length_Vz)

% Naiquist theorem (Kotelnikov)
t_wats=linspace(1, Length_Vz, L);
td_wats=linspace(1, Length_Vz, Length_Vz);
s_wats=Vz.';
d_wats=[td_wats' s_wats'];
y_wats=pulstran(t_wats, d_wats, 'sinc'); % Vz interpolated in all WATS data

    Vz_interpolated=y_wats(1:length(y_wats)); % change work boundaries of vector for each case
end