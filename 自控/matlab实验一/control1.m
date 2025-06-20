s=tf("s");
G1=1/(s+1);
G2=1/(s+2);
G3=1/(s+3);
G4=-1/(s+4);
G5=1/(s+5);
G6=-1;
sys=append(G1,G2,G3,G4,G5,G6);
Q=[1 5 6;
    2 1 4;
    3 2 0;
    4 3 0;
    5 2 0;
    6 3 0];
inputs=1;
outputs=3;
sysc=connect(sys,Q,inputs,outputs);
G=tf(sysc)