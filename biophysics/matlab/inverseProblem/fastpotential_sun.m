function V = fastpotential_sun(r, m, s)

% /**********************************************************************/
% /* Sun M., "An Efficient Algorithm for Computing Multishell Spherical */
% /* Volume Conductor Models in EEG Dipole Source Localization",        */
% /* IEEE Trans. Biomed. Eng., 44(12), pp. 1243--52, 1997.              */
% /*                                                                    */
% /* Code copied directly from page 1248 of above publication!          */
% /*                                                                    */
% /* NB! Code expects and returns the following units:                  */
% /*     Dipole location    (r): cm                                     */
% /*     Dipole moment      (m): mA cm                                  */
% /*     Electrode location (s): cm                                     */
% /*     *** Output voltage    : mV                                     */
% /*                                                                    */
% /* The (default) 4-shell conductor model is by Stok (86):             */
% /*                                                                    */
% /*               Conductivity (S/m)     Radius (relative to scalp, R) */
% /*     Brain          0.33                0.8400  (63 mm)             */
% /*     CSF            1.0                 0.8667  (65 mm)             */
% /*     Skull          0.0042              0.9467  (71 mm)             */
% /*     Scalp          0.33                1.0000  (75 mm)             */
% /*                                                                    */
% /*                                                                    */
% /* C. Bailey, 17/06/02                                                */
% /**********************************************************************/
% double FastPotential (double *r, double *m, double *s) {
% 
%   int i;
%   double R, R2,
%     Q, Q2, Q3, Q5, Q7,
%     x, x2, f, f2, f12, fx,
%     vf, vr, p, q, t0, r0, T[3],

C = 24.114385316954;

%     /* These are for model of Cuffin & Cohen -79 */
%cn = [-3.45393858774666, -0.81750157854542, -0.09602868634859, 0.12927294546578];
%a  = [3.71904029916523, -0.27826260492429, 0.01336227893567, -0.00020138542995];
%     
%     /* These are for model of Stok -86 (PhD Thesis) */
cn = [0.0, -0.06141320805, 0.05342974625, 0.02582167642];
a  = [2.094413085618, -0.212973206782, 0.015140880154, -0.000311147671];

R2 = dot(s,s);
R = sqrt(R2);
p = dot(r,r);
f = sqrt(p);
q = dot(r,s);

if (f == 0.0),
    V = 44.2464 *dot(m,s) / (R2*R);
    return
end

% for (i = 0; i < 3; i++)
%     T[i] = s[i]*p - r[i]*q;
T = s*p - r*q;

x = dot(T,T);

if (x == 0.0)
    t0 = 0;
else
    t0 = dot(m,T) / sqrt(x);
end

r0 = dot(m,r) / f;
x = q / (f*R);

foo = f/R;
f = foo;

f2 = f*f;
x2 = x*x;
fx = f*x;
f12 = 1.0-f2;

Q2 = 1.0/(1.0 - 2.0*fx + f2);
Q = sqrt(Q2);
Q3 = Q2*Q;
Q5 = Q3*Q2;
Q7 = Q5*Q2;

vr = cn(2)*x + cn(3)*f*(1.5*x2 - 0.5) + a(1)*(Q - 1.0)/f + a(2)*(x-f)*Q3 + ...
    a(3)*(f*(x2 + f2 - 2.0) + x*f12)*Q5 + cn(4)*f2*x*(2.5*x2 - 1.5) + ...
    a(4)*(x + f*((5.0*x2 - 4.0) + f*((x*(x2 - 9.0)) + f*((10.0 - 2.0*x2)-fx-f2))))*Q7;

vf = sqrt(1-x2)*(cn(2) + 3.0*fx*cn(3)/2.0 + a(1)*Q*(1.0 + Q)/(Q - fx*Q + 1.0) + ...
    a(2)*Q3 + a(3)*Q5*(f12-f2+fx) + cn(4)*(5.0*x2 - 1.0)*f2/2.0 + ...
    a(4)*(1.0 + f2*(4.0*f2 - fx + x2 - 10.0) + 5.0*fx)*Q7);


V = C*(t0*vf + r0*vr)/R2;

    
