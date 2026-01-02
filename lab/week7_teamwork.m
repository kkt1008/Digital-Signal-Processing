% 1번 
B=9;
% a) 0.12345
a=0.12345; 
a_s=a*2^B;
ab = dec2bin(fix(a_s),B+1);

B=9;
% b) -0.56789
b = -0.56789;
b_s=b*2^B;
% 부호 비트를 맨 뒤에 비트에 1을 추가하는 방식이라 출력 결과가 계산과 다름
% 부호-크기 형식으로 출력하기 위해 살짝의 조작 필요
% 부호-크기 표현으로 10비트 이진수 생성
%bin_abs = dec2bin(abs_val, B); % 9비트 양수 이진수: '100100010'
%sign_bit = (val < 0); % 음수면 1, 양수면 0
%bin_sign_mag = [num2str(sign_bit), bin_abs]; % 부호 비트 + 9비트

bb = dec2bin(fix(b_s),B); % dec2bin 함수는 음수일 경우 16bit의 2의보수로 변환
% 16bit 중 뒤에서 10개 bit만 추출
%bb = bb(end-:end); % 9비트 양수 이진수: '100100010'
disp(bb);
sign_bit = (b_s < 0); % 음수면 1, 양수면 0
bb = [num2str(sign_bit), bb]; % 부호 비트 + 9비트
disp(bb); % 정답: 1 100100010


% c) 0.38452386
c = 0.38452386;
c_s = c*2^B;
cb = dec2bin(fix(c_s),B+1);

% d) -0.762349
d = -0.762349;
d_s = d*2^B;
db = dec2bin(fix(d_s),B);
db = db(end-9:end);
sign_bit = (b_s < 0); % 음수면 1, 양수면 0
db = [num2str(sign_bit), db]; % 부호 비트 + 9비트


% e) -0.90625
e = -0.90625;
e_s = e*2^B;
eb = dec2bin(fix(e_s),B);
eb = eb(end-9:end);
sign_bit = (b_s < 0); % 음수면 1, 양수면 0
bb = [num2str(sign_bit), eb]; % 부호 비트 + 9비트


disp('10-bit sign-magnitude');
disp(['a 10-bit sign: ',ab]);
disp(['b 10-bit sign: ',bb]);
disp(['c 10-bit sign: ',cb]);
disp(['d 10-bit sign: ',db]);
disp(['e 10-bit sign: ',eb]);

% a-1의 보수변형
%a1 = OnesComplement(ab,B+1);
% b-1의 보수변형
b1 = OnesComplement(-bb,B+1);
% c-1의 보수변형

% d-1의 보수변형

% e-1의 보수변형

