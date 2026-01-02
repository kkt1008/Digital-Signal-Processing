B=9;
% a) 0.12345
a = 0.12345​
a_s = a*2^B; % approxi 63 0 ​
ab = dec2bin(fix(a_s),B);​

% b) -0.56789
b = 0.56789;​
b_s = b*2^B; % approxi 290​
bb = dec2bin(fix(b_s),B); % 1 ​

% c) 0.38452386
c = 0.38452386;​
c_s = c*2^B; % approxi 290​
cb = dec2bin(fix(c_s),B); % 0 ​

% d) 0762349
d = 0.762349;​
d_s = d*2^B; % approxi 290​
db = dec2bin(fix(d_s),B); % 1 ​

% e) 0.90625​
e = 0.90625;​
e_s = e*2^B; % approxi 290​
eb = dec2bin(fix(e_s),B); % 1 ​

disp(['a 10-bit sign : 0.',ab]);​
disp(['b 10-bit sign : 1.',bb]);​
disp(['c 10-bit sign : 0.',cb]);​
disp(['d 10-bit sign : 1.',db]);​
disp(['e 10-bit sign : 1.',eb]);​

​

% a1 = OnesComplement(-bb,B+1);​
b2 = TwosComplement(-fix(b_s),B+1);​
b2s = dec2bin(b2,B);​

d2 = TwosComplement(-fix(d_s),B+1);​
d2s = dec2bin(d2,B);​

e2 = TwosComplement(-fix(e_s),B+1);​
e2s = dec2bin(e2,B);​

disp(['b 2s comp : 1.',b2s]);​
disp(['d 2s comp : 1.',d2s]);​
disp(['e 2s comp : 1.',e2s]);​