#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White transmit 1.0}
camera {orthographic
  right -25.94*x up 12.91*y
  direction 1.00*z
  location <0,0,50.00> look_at <0,0,0>}
light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}

#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.070;
#declare Rbond = 0.100;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
      translate LOC}
#end

cylinder {< -5.33, -12.81,  -3.46>, <  3.07, -12.49,  -2.92>, Rcell pigment {Black}}
cylinder {< -0.70, -11.66, -10.40>, <  7.70, -11.33,  -9.86>, Rcell pigment {Black}}
cylinder {< -1.73,  10.00,  -7.48>, <  6.67,  10.32,  -6.94>, Rcell pigment {Black}}
cylinder {< -6.35,   8.84,  -0.54>, <  2.04,   9.17,   0.00>, Rcell pigment {Black}}
cylinder {< -5.33, -12.81,  -3.46>, < -0.70, -11.66, -10.40>, Rcell pigment {Black}}
cylinder {<  3.07, -12.49,  -2.92>, <  7.70, -11.33,  -9.86>, Rcell pigment {Black}}
cylinder {<  2.04,   9.17,   0.00>, <  6.67,  10.32,  -6.94>, Rcell pigment {Black}}
cylinder {< -6.35,   8.84,  -0.54>, < -1.73,  10.00,  -7.48>, Rcell pigment {Black}}
cylinder {< -5.33, -12.81,  -3.46>, < -6.35,   8.84,  -0.54>, Rcell pigment {Black}}
cylinder {<  3.07, -12.49,  -2.92>, <  2.04,   9.17,   0.00>, Rcell pigment {Black}}
cylinder {<  7.70, -11.33,  -9.86>, <  6.67,  10.32,  -6.94>, Rcell pigment {Black}}
cylinder {< -0.70, -11.66, -10.40>, < -1.73,  10.00,  -7.48>, Rcell pigment {Black}}
atom(< -5.68,  -5.39,  -2.46>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #0 
atom(< -2.88,  -5.28,  -2.28>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #1 
atom(< -0.08,  -5.17,  -2.10>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #2 
atom(< -4.14,  -5.00,  -4.77>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #3 
atom(< -1.34,  -4.89,  -4.59>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #4 
atom(<  1.46,  -4.79,  -4.41>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #5 
atom(< -2.59,  -4.62,  -7.09>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #6 
atom(<  0.20,  -4.51,  -6.90>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #7 
atom(<  3.00,  -4.40,  -6.72>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #8 
atom(< -4.34,  -2.95,  -2.86>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #9 
atom(< -1.54,  -2.85,  -2.68>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #10 
atom(<  1.26,  -2.74,  -2.50>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #11 
atom(< -2.80,  -2.57,  -5.18>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #12 
atom(<  0.00,  -2.46,  -5.00>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #13 
atom(<  2.80,  -2.35,  -4.81>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #14 
atom(< -1.25,  -2.18,  -7.49>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #15 
atom(<  1.54,  -2.08,  -7.31>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #16 
atom(<  4.34,  -1.97,  -7.13>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #17 
atom(<  2.61,  -0.34,  -2.91>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #18 
atom(< -2.98,  -0.53,  -3.28>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #19 
atom(< -0.20,  -0.41,  -3.08>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #20 
atom(<  4.16,   0.06,  -5.23>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #21 
atom(< -1.46,  -0.12,  -5.55>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #22 
atom(<  1.35,  -0.06,  -5.41>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #23 
atom(<  5.66,   0.49,  -7.54>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #24 
atom(<  0.08,   0.19,  -7.89>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #25 
atom(<  2.91,   0.35,  -7.72>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #26 
atom(< -1.36,   2.56,  -8.48>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #27 
atom(<  1.38,   2.60,  -8.28>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #28 
atom(< -0.39,   1.81,  -1.19>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #29 
atom(< -4.42,   1.79,  -3.86>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #30 
atom(< -1.62,   2.09,  -3.63>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #31 
atom(<  1.23,   1.97,  -3.53>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #32 
atom(< -2.89,   2.19,  -6.17>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #33 
atom(< -0.10,   2.28,  -5.99>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #34 
atom(<  2.72,   2.39,  -5.82>, 1.21, rgb <0.82, 0.82, 0.87>, 0.0, ase2) // #35 
atom(< -1.18,   3.46,  -2.16>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #36 
atom(< -0.29,   4.63,  -2.46>, 0.68, rgb <0.56, 0.56, 0.56>, 0.0, ase2) // #37 
atom(< -2.11,   3.79,  -1.66>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #38 
atom(< -0.82,   5.33,  -3.13>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #39 
atom(<  0.67,   4.36,  -2.94>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #40 
atom(< -0.06,   5.18,  -1.53>, 0.28, rgb <1.00, 1.00, 1.00>, 0.0, ase2) // #41 
