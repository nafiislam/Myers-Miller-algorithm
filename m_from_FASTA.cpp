#include <iostream>
#include <string>
#include <vector>
#include <climits>
#include <algorithm>
#include <fstream>
using namespace std;

int MATCH, MISMATCH; //replacement costs = 0, 1
float G, H; // for gap cost = G + H*k, G=2, H=0.5

string A, B; // input strings
vector<int> CC, DD, RR, SS; // 1D arrays for cost calculation
vector<int> S; // output alignment-mapper
int Slen; // number of entries in S

string readFasta(string filename){
    ifstream fin(filename);
    string line, seq = "";

    while(getline(fin, line)){
        if(line.empty() || line[0] == '>'){
            continue;
        }
        if(line.back() == '\r'){ // Windows line ending
            line.pop_back();
        }
        seq = seq + line;
    }
    return seq;
}

// Cost function
int w(char a, char b){
    if(a == b){
        return MATCH;
    }
    else{
        return MISMATCH;
    }
}

int minimum(int x, int y){
    if(x < y){
        return x;
    }
    else{
        return y;
    }
    
}

int gapcost(int k){
    if(k == 0){
        return 0;
    }
    else{
        return G + H * k;
    }
}

// 0-> match or replace, -something -> delete in B, +something -> insert in B
void addop(int op){
    if(op != 0 && Slen > 0){ //checking with the prev values if it is gap as well
        int last = S[Slen - 1];
        if((last > 0 && op > 0) || (last < 0 && op < 0)){
            S[Slen - 1] = last + op;
            return;
        }
    }
    S[Slen] = op;
    Slen = Slen + 1;
}

// Cost-only linear-space sweep (Fig 1B)
// a0, b0 are offsets into A, B; M, N are lengths of the substrings
// reverse=0 -> strings are read left to right; otherwise right to left
// tb = gap-initiation cost (0 if bracketed by an adjacent delete, G otherwise)
void sweep(int a0, int b0, int M, int N, int rows, int tb, int reverse,
           vector<int> &C, vector<int> &D){
    int i, j;
    int t, s, c, e;
    char ai, bj;

    C[0] = 0;
    t = G;
    for(j = 1; j <= N; j++){
        t = t + H;
        C[j] = t;
        D[j] = t + G;
    }
    t = tb; // the starred line [*], later fixed
    for(i = 1; i <= rows; i++){
        if(reverse == 0){
            ai = A[a0 + i - 1];
        } 
        else{
            ai = A[a0 + M - i];
        }              
        s = C[0];
        t = t + H;
        c = t;
        C[0] = c;
        e = t + G;
        for(j = 1; j <= N; j++){
            if(reverse == 0){
                bj = B[b0 + j - 1];
            }
            else{
                bj = B[b0 + N - j];
            }
            e = minimum(e, c + G) + H;
            D[j] = minimum(D[j], C[j] + G) + H;
            c = minimum(D[j], e);
            c = minimum(c, s + w(ai, bj));
            s = C[j];
            C[j] = c;
        }
    }
    D[0] = C[0]; // from paper corrections
}

// tb / te = gap-initiation cost delete at the start / end conversion (G normally, 0 when adjacent delete each).
void diff(int a0, int b0, int M, int N, int tb, int te){
    int j, mid, jstar, type1;
    int best, c1, c2, c, mt;

    if(N == 0){
        if(M > 0){
            addop(-M); // delete in B (gap in B)
        }
    }
    else if(M == 0){
        addop(N); // insert in B
    }
    else if(M == 1){
        // minimum cost calculation
        mt = minimum(tb, te);
        best = mt + H + gapcost(N); // delete a1 (gap in B), insert B
        jstar = 0;
        for(j = 1; j <= N; j++){
            c = gapcost(j - 1) + w(A[a0], B[b0 + j - 1]) + gapcost(N - j);
            if(c <= best){
                best = c;
                jstar = j;
            }
        }

        // writing it into S
        if(jstar == 0){
            if(tb < te){ // make gaps of B to the left
                addop(-1);
                addop(N);
            }
            else{ // make gaps of B to the right
                addop(N);
                addop(-1);
            }
        } 
        else{
            if(jstar - 1 > 0){ // safety check
                addop(jstar - 1);
            }
            addop(0); // match or replace
            if(N - jstar > 0){ // safety check
                addop(N - jstar);
            }
        }
    }
    else{
        mid = M / 2; // middle point of A

        sweep(a0, b0, M, N, mid, tb, 0, CC, DD); // forward sweep
        sweep(a0, b0, M, N, M - mid, te, 1, RR, SS); // reverse sweep

        // find the optimal midpoint in B
        best = INT_MAX;
        jstar = 0;
        type1 = 1;
        for(j = 0; j <= N; j++){
            c1 = CC[j] + RR[N - j];
            c2 = DD[j] + SS[N - j] - G;
            if(c1 < best){
                best = c1;
                jstar = j;
                type1 = 1;
            }
            if(c2 < best){
                best = c2;
                jstar = j;
                type1 = 0;
            }
        }

        if(type1 == 1){
            diff(a0, b0, mid, jstar, tb, G);
            diff(a0 + mid, b0 + jstar, M - mid, N - jstar, G, te);
        } else{
            diff(a0, b0, mid - 1, jstar, tb, 0);
            addop(-2); // delete a_i* and a_i*+1 (create 2 gaps in B)
            diff(a0 + mid + 1, b0 + jstar, M - mid - 1, N - jstar, 0, te);
        }

    }
}

// Cost of the alignment
int alignmentCost(int len){
    int k, i = 0, j = 0, op;
    int total = 0;
    for(k = 0; k < len; k++){
        op = S[k];
        if(op == 0){ // match or replace
            total = total + w(A[i], B[j]);
            i++;
            j++;
        }
        else if(op < 0){ // gaps in B
            total = total + gapcost(-op);
            i = i - op; // going right of A
        }
        else{ // insert in B
            total = total + gapcost(op);
            j = j + op; // going right of B
        }
    }
    return total;
}

// Compute the minimum cost and fill S / Slen.
int align(){
    int M = A.size();
    int N = B.size();
    CC.assign(N + 1, 0);
    DD.assign(N + 1, 0);
    RR.assign(N + 1, 0);
    SS.assign(N + 1, 0);
    S.assign(M + N + 2, 0); // maximum possible length of S is M + N, +2 for safety
    Slen = 0;
    diff(0, 0, M, N, G, G);
    return alignmentCost(Slen);
}

void display(ostream &out){
    string top = "", mid = "", bot = "";
    int k, t, i = 0, j = 0, op;
    for(k = 0; k < Slen; k++){
        op = S[k];
        if(op == 0){ // match or replace
            top = top + A[i];
            bot = bot + B[j];
            if(A[i] == B[j]){ // match
                mid = mid + '!';
            }
            else{ // replace
                mid = mid + '|';
            }
            i++;
            j++;
        }
        else if(op < 0){// gaps in B
            for(t = 0; t < -op; t++){
                top = top + A[i];
                mid = mid + ' ';
                bot = bot + '-';
                i++;
            }
        } else{ // insert in B
            for(t = 0; t < op; t++){
                top = top + '-';
                mid = mid + ' ';
                bot = bot + B[j];
                j++;
            }
        }
    }
    out << top << endl << mid << endl << bot << endl;
}

int main(int argc, char *argv[]){
    // Example 1 from the paper (mismatch 1, match 0, gap(k) = 2 + 0.5k),
    // But as we are using int, we need to double. So, cost will be 8 instead of 4.
    MATCH = 0;
    MISMATCH = 2; //1
    G = 4; //2
    H = 1; //0.5
 
    A = readFasta(argv[1]);
    B = readFasta(argv[2]);


    int cost = align();
    // cout << "cost is " << cost << endl;
    cout << "Alignment-mapper:";
    for(int k = 0; k < Slen; k++){
        cout << " " << S[k];
    }
    cout << endl;

    ofstream fout(argv[3]);

    fout << "Length A = " << A.size() << ", length B = " << B.size() << endl;
    fout << "Alignment cost = " << cost << endl << endl;
    display(fout);
    fout.close();

    return 0;
}