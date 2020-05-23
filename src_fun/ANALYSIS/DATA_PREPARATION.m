%prepare sim data
weights = h5read('MODEL_LARGE.hdf5','/weights');
STATE_TEMP = h5read('MODEL_LARGE.hdf5','/STATE_TEMP');

S1 = size(STATE_TEMP,1);
S2 = size(STATE_TEMP,2);
S3 = size(STATE_TEMP,3);

ST=[];

for i = 1:S3
    ST = [ST;transpose(STATE_TEMP(:,:,i))];
end

ENS = zeros(S2*S3,S3);

ENS(1:6,1)=1;

for i = 2:S3
    A = (((i-1)*6)+1);
    B = (i*6);
    ENS(A:B,i)=1;
end

ST_ENS = [ST ENS];

SHUF_ENS = ST_ENS(randperm(size(ST_ENS,1)),:);

data = [SHUF_ENS; SHUF_ENS; SHUF_ENS];
UDF = data(1:500,(S1+1):end);
data = data(1:500,1:S1);

coords = randi(512,S1,2);
coords = [coords randi(10,S1,1)];

data(data<1)=0;

filename = 'MODEL_LARGE.mat';
clear A B ENS i S1 S2 S3 SHUF_ENS ST ST_ENS
save(filename);