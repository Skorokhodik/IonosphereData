% About method of Gravity Wave selection

% Here we work with two instruments from Dynamics Explorer 2 satellite
- NACS and WATS 

NACS
% First step -	4-7 
determine date what you want to research
year with 364/365 day

% Second step -	9-22
download NACS data and determine parameters: 
UT, Satellite Altitude, Latitude, Longitude, LST, Oxigen cons.,
N2 cons., He cons., Ar cons., Orbit number

% Third step -	25-33
looking for breaking point in NACS
- msec to sec
- Looking for break points and put in massive [UT(k-1) UT(k) BreakeSeconds]
%BreakeSeconds - how many seconds left between UT(k-1) and UT(k)

% Fourth step -	36-53
Oxigen concentration massive 
use function GravitationWave for detect it in Oxigen data set
				55-64
Nitrogen concentration massive 
use function GravitationWave for detect it in Nitrogen data set

% Fifth step -	66-95
Calculate NOISE 
	1 - spectra NOISE
	2 - iFFT of NOISE
	3 - Histogram of NOISE

WATS
% Sixth step - 97-
		1 - way to WATS data
		2 - download WATS data of the same date what NACS



%% PICTURES
44-65
Figure1,	1 - Ox (red line) and N2(blue line) concentration 
			2 - GW in Ox (red line) and N2(blue line)
77-93
Figure2,	1 - NOISE spectra
			2 - NOISE in space
			3 - histogram of NOISE

%% FUNCTIONs
	GravitationWave
		1 - take parameter massive
		2 - build spectra
		3 - select GW in spectra
		
	WATS_data_transmission
		1 - 
		2 - 
		
	Vertical_wind_WATS
		1 - 
		2 - 
		
	Wind_Interpolated_WATS
		Looking for aria where WATS & NACS data agreement
		and interpolated wind data from WATS
		1 - 
		2 - 
		3 - 
		
		Naiquist_theorem
		
		