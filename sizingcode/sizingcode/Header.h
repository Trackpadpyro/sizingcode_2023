#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h> 

using namespace std;

float lambda = 0.7; // battery max discharge
float t_length = 1000; // track length, regarding this, I couldn't find a suitable estimate for turn radius 
//(most formulae online just lead to nonsense and impossible configurations), 
//so I looked at a couple past comps and noticed that approximately for every second in a turn, 
//most planes spend 2s in the straight, and hence I used 1.5x609.6 ~ 900m, and rounded up to 1000m for the estimate for track length.
float passenger = 0.038464632; // passenger weight

float m2_max = 0.077052;
float m3_max = 0.0214029;

int constraint_num = 9;
int variable_num = 9;
vector<float> constraint_list = { 12.5, 25, 0.68, 0.68, 2.322576, 1 / (120 * lambda), 6.096, 6.096, 0 };
vector<float> u_bound = { 10        ,2.32       ,25         ,25         ,3.8464632  ,2800       ,2800        ,1.6        ,5 };
vector<float> l_bound = { 1         ,0.4        ,0.6        ,0          ,0          ,0.33       ,0.33        ,0.001      ,0.002 };
// 0 => constraint list < contraint; 1 => constraint list > contraint
vector<bool> constraint_sgn = { 0, 1, 1, 1, 1, 0, 1, 1, 0 };

float magsq(vector<float> v);
vector<float> add(vector<float> a, vector<float> b);
vector<float> scalar(vector<float> a, float b);
vector<float> normalize(vector<float> v, float step_size);
float constraint(vector<float> o, int n);
bool check_valid(vector<float> a, int n);
vector<float> linearization(vector<float> v, int n, float w);
vector<float> perp(vector<float> a, vector<float> c);
bool check_valid_all(vector<float> x);

float magsq(vector<float> v) {
    float sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += (float)pow(v[i], 2);
    }
    return sum;
}

vector<float> add(vector<float> a, vector<float> b) {
    if (a.size() != variable_num || b.size() != variable_num) {
        cout << "invalid vector (1)" << endl;
        return { -1 };
    }
    else {
        vector<float> temp;
        for (int i = 0; i < variable_num; i++) {
            temp.push_back(a[i] + b[i]);
        }
        return temp;
    }
}

vector<float> scalar(vector<float> a, float b) {
    if (a.size() != variable_num) {
        cout << "invalid vector (2)" << endl;
        return { -1 };
    }
    else {
        vector<float> temp;
        for (int i = 0; i < variable_num; i++) {
            temp.push_back(a[i] * b);
        }
        return temp;
    }
}

vector<float> normalize(vector<float> v, float step_size) {
    vector<float> result;
    for (int i = 0; i < v.size(); i++) {
        result.push_back((step_size * v[i]) / sqrt(magsq(v)));
    }
    return result;
}

float constraint(vector<float> o, int n) {
    float ka, kt;
    if (o.size() != variable_num) {
        cout << "invalid vector (3)" << endl;
        return -1;
    }
    else {
        switch (n) {
        case 1: // must complete M1 in time
            return sqrt(o[2] / (o[1] * o[7] * 0.613));
        case 2: // total weight < 55lb
            return o[2] + o[3];
        case 3: // thrust exceeds drag for M2
            return sqrt((o[2] + o[3]) / (o[1] * o[7] * 0.613)) * (o[2] + o[3]) * o[8] / o[7] / o[5];
        case 4: // thrust exceeds drag for M3
            return sqrt((o[2] + o[4]) / (o[1] * o[7] * 0.613)) * (o[2] + o[4]) * o[8] / o[7] / o[6];
        case 5: // wing area
            return o[0] * o[1];
        case 6: // range is at least 3 laps
            return sqrt((o[2] + o[3]) / (o[1] * o[7] * 0.613)) / o[5];
        case 7: // takeoff distance for M2
            ka = o[1] * (0.04 * 1.5 - 0.1) / (0.612 * (o[2] + o[3]));
            kt = o[5] * 0.68 / (sqrt((o[2] + o[3]) / (o[1] * 1.5 * 0.613)) * (o[2] + o[3])) - 0.04;
            return log(1 + ka * 1.44 * (o[2] + o[3]) / (o[1] * 1.5 * 0.613 * kt)) / (2 * 9.18 * ka);
        case 8: // takeoff distance for M3
            ka = o[1] * (0.04 * 1.5 - 0.1) / (0.612 * (o[2] + o[4]));
            kt = o[6] * 0.68 / (sqrt((o[2] + o[4]) / (o[1] * 1.5 * 0.613)) * (o[2] + o[4])) - 0.04;
            return log(1 + ka * 1.44 * (o[2] + o[4]) / (o[1] * 1.5 * 0.613 * kt)) / (2 * 9.18 * ka);
        case 9: // min cd cl contraint
            return ((float)pow(o[7], 2)) / (o[0] * 3.1415926 * 0.7) - o[8];
        }
    }
}

vector<float> linearization(vector<float> v, int n, float w) {
    vector<float> approx;
    vector<float> temp;

    if (v.size() != variable_num) {
        cout << "invalid vector (4)" << endl;
        return { -1 };
    }
    for (int i = 0; i < variable_num; i++) {
        temp = v;
        temp[i] += w;
        approx.push_back((constraint(temp, n) - constraint(v, n)) / w);
    }
    return approx;
}

vector<float> perp(vector<float> a, vector<float> c) {
    if (a.size() != variable_num || c.size() != variable_num) {
        cout << "invalid vector (5)" << endl;
        return { -1 };
    }
    else {
        vector<float> b = normalize(c, 1);
        float temp = 0;
        for (int i = 0; i < variable_num; i++) {
            temp -= a[i] * b[i];
        }
        return add(a, scalar(b, temp));
    }
}

bool check_valid(vector<float> a, int n) {
    if (a.size() != variable_num || n >= constraint_num || n < 0) {
        cout << "invalid vector (6)" << endl;
        return false;
    }
    else if (constraint(a, n + 1) <= constraint_list[n] && constraint_sgn[n]) {
        //cout << constraint(a, n + 1) << " leq " << constraint_list[n] << endl;
        return true;
    }
    else if (constraint(a, n + 1) >= constraint_list[n] && !constraint_sgn[n]) {
        //cout << constraint(a, n + 1) << " geq " << constraint_list[n] << endl;
        return true;
    }
    else {
        return false;
    }
}

bool check_valid_all(vector<float> x) {
    bool valid = 1;
    if (x.size() != variable_num) {
        cout << "invalid vector (10)" << endl;
        return 0;
    }
    else {
        for (int i = 0; i < x.size(); i++) {
            if (x[i] > u_bound[i] || x[i] < l_bound[i]) {
                valid = 0;
            }
        }
        for (int i = 0; i < constraint_num; i++) {
            if (!check_valid(x, i)) {
                valid = 0;
                break;
                //cout << i << ": " << constraint(x, i) << endl;
            }
        }
    }
    return valid;
}

class Parameter_Set {
public:
    float asp_ratio;    //0
    float s_area;       //1
    float weight;       //2
    float p2_weight;    //3
    float p3_weight;    //4
    float power2;       //5
    float power3;       //6
    float cl;           //7
    float cd;           //8

    Parameter_Set(float x_1, float x_2, float x_3, float x_4, float x_5, float x_6, float x_7, float x_8, float x_9) {
        asp_ratio = x_1;
        s_area = x_2;
        weight = x_3;
        p2_weight = x_4;
        p3_weight = x_5;
        power2 = x_6;
        power3 = x_7;
        cl = x_8;
        cd = x_9;
    }
    Parameter_Set(vector<float> x) {
        if (x.size() != variable_num) {
            cout << "invalid vector (7)" << endl;
        }
        else {
            asp_ratio = x[0];
            s_area = x[1];
            weight = x[2];
            p2_weight = x[3];
            p3_weight = x[4];
            power2 = x[5];
            power3 = x[6];
            cl = x[7];
            cd = x[8];
        }
    }
    float score_M2() {
        return p2_weight * vel(p2_weight) / (3 * t_length);
    }
    float score_M3() {
        return (lambda * p3_weight * vel(p3_weight)) / (t_length * passenger * power3);
    }
    float score() {
        return score_M2() / m2_max + score_M3() / m3_max + 4;
    }
    float vel(float w) {
        return sqrt((weight + w) / (s_area * cl * 0.613));
    }
    vector<float> gradient_M2() {
        vector<float> gradient_M2 =
        { 0,
            score_M2() * ((float)-0.5) / s_area,
            score_M2() / (2 * (weight + p2_weight)),
            vel(p2_weight) / (3 * t_length) * (((float)1) + ((float)0.5) * (p2_weight / (p2_weight + weight))),
            0,
            0,
            0,
            score_M2() * ((float)-0.5) / cl,
            0
        };
        return gradient_M2;
    }
    vector<float> gradient_M3() {
        vector<float> gradient_M3 =
        { 0,
            score_M3() * ((float)-0.5) / s_area,
            score_M3() / (2 * (weight + p3_weight)),
            0,
            ((lambda * vel(p3_weight)) / (t_length * passenger * power3)) * (1 + ((float)0.5) * (p3_weight / (p3_weight + weight))),
            0,
            score_M3() * ((float)-1) / (power3),
            score_M3() * ((float)-0.5) / cl,
            0
        };
        return gradient_M3;
    }
    vector<float> step(float step_size, int state) {
        vector<float> current = { asp_ratio, s_area, weight, p2_weight, p3_weight, power2, power3, cl, cd };
        
        vector<float> gradient;
        switch (state) {
        case 1:
            gradient = add(scalar(gradient_M3(), (float)(1 / m3_max)), scalar(gradient_M2(), (float)(1 / m2_max)));
            break;
        case 2:
            gradient = gradient_M2();
            break;
        case 3:
            gradient = gradient_M3();
            break;
        }

        vector<float> s = scalar(gradient, step_size);

        for (int i = 0; i < s.size(); i++) {
            if (s[i] + current[i] > u_bound[i]) {
                s[i] = u_bound[i] - current[i];
            }
            if (s[i] + current[i] < l_bound[i]) {
                s[i] = l_bound[i] - current[i];
            }
        }
        return s;
    }
    void paste(int state) {
        cout << "asp_ratio: " << asp_ratio << "; "
            << "s_area: " << s_area << "; "
            << "weight: " << weight << "; "
            << "p2_weight: " << p2_weight << "; "
            << "p3_weight: " << p3_weight << "; "
            << "power2: " << power2 << "; "
            << "power3: " << power3 << "; "
            << "cl: " << cl << "; "
            << "cd: " << cd << endl;
        switch (state) {
        case 1:
            cout << "score: " << score() << "; score M2: " << score_M2() << "; score M3: " << score_M3() << endl;
            break;
        case 2:
            cout << "score M2: " << score_M2() << endl;
            break;
        case 3:
            cout << "score M3: " << score_M3() << endl;
            break;
        }
    }
    vector<float> status() {
        vector<float> current = { asp_ratio, s_area, weight, p2_weight, p3_weight, power2, power3, cl, cd };
        return current;
    }
    void update(vector<float> s) {
        if (s.size() == variable_num) {
            asp_ratio += s[0];
            s_area += s[1];
            weight += s[2];
            p2_weight += s[3];
            p3_weight += s[4];
            power2 += s[5];
            power3 += s[6];
            cl += s[7];
            cd += s[8];
        }
        else {
            cout << "invalid vector (8)" << endl;
        }
    }
    void set(vector<float> s) {
        if (s.size() == variable_num) {
            asp_ratio = s[0];
            s_area = s[1];
            weight = s[2];
            p2_weight = s[3];
            p3_weight = s[4];
            power2 = s[5];
            power3 = s[6];
            cl = s[7];
            cd = s[8];
        }
        else {
            cout << "invalid vector (9)" << endl;
        }
    }
    void set(Parameter_Set s) {
        asp_ratio = s.asp_ratio;
        s_area = s.s_area;
        weight = s.weight;
        p2_weight = s.p2_weight;
        p3_weight = s.p3_weight;
        power2 = s.power2;
        power3 = s.power3;
        cl = s.cl;
        cd = s.cd;
    }

    void randomize(int part_num) {
        vector<float> f;
        while (true) {
            f.clear();
            for (int i = 0; i < variable_num; i++) {
                f.push_back((u_bound[i] - l_bound[i]) * ((float)(rand() % part_num)) / ((float)part_num) + l_bound[i]);
            }
            if (check_valid_all(f)) {
                set(f);
                //cout << 1 << endl;
                break;
            }
        }
    }
};
