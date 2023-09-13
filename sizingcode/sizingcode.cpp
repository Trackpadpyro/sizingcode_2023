#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h> 
#include <fstream>
#include <string>

#include "Header.h"

using namespace std;

int part_num = 100000;
int runtime = 200000000;

Parameter_Set descent(Parameter_Set t1, float max_step_size, float factor, int n, int state, bool silent) {
    vector<float> temp1;
    float temp3 = 0;
    float score;

    if (!silent) { 
        cout << "Start: " << endl; 
        t1.paste(state);
    }
    for (int j = 0; j < n; j++) {
        temp1 = t1.step(max_step_size, state);
        for (int i = 0; i < constraint_num; i++) {
            if (!check_valid(add(t1.status(), temp1), i)) {
                temp1 = perp(temp1, linearization(t1.status(), i, max_step_size / 100));
            }
        }
        while (true) {
            if (check_valid_all(add(temp1, t1.status()))) {
                break;
            }
            else {
                temp1 = scalar(temp1, factor);
            }
        }
        t1.update(temp1);
        score = t1.score(state);
        if (j % 10 == 0 && j != 0) {
            if (temp3 == score) {
                break;
            }
            temp3 = score;
        }
    }
    if (!silent) {
        cout << "End: " << endl;
        t1.paste(state);
    }
    return t1;
}
int main()
{
    srand(time(NULL));
    vector<float> c = { 3.15096, 0.638872, 44.4805, 21.8988, 11.9057, 751.007, 587.023, 0.114193, 0.0218818 };
    

    Parameter_Set x(c);
    Parameter_Set max(c);
    int mission = 1;
    max.paste(1);
    max.paste(10);

    ofstream rec1("info_M1.txt");
    for (int i = 1; i <= runtime; i++) {
        x.randomize(part_num);
        if (max.score(mission) < x.score(mission)) {
            max.set(descent(x, 0.5, 0.5, 10000, mission, 1));
            for (int i = 0; i < x.status().size(); i++) {
                rec1 << max.status()[i] << ", ";
            }
            rec1 << max.score() << ", " << max.score_M2() << ", " << max.score_M3() << endl;
        }

        for (int i = 0; i < x.status().size(); i++) {
            rec1 << x.status()[i] << ", ";
        }
        rec1 << x.score() << ", " << x.score_M2() << ", " << x.score_M3() << endl;

        if (i % 1000 == 0) {
            max.paste(1);
        }
    }
    rec1.close();
    max.paste(1);
    max.paste(10);
    exit(runtime);
}
/*
asp_ratio: 3.63718; s_area: 0.406715; weight: 43.0023; p2_weight: 73.899; p3_weight: 12.2799; power2: 897.36; power3: 224.06; cl: 0.349859; cd: 0.0353029
score: 5.72063; score M2: 0.90179; score M3: 9.21449
M1 cruise speed: 22.2036; M2 cruise speed: 36.609; M3 cruise speed: 25.1751
M1 cruise drag: 4.33919; M2 cruise drag: 11.796; M3 cruise drag: 5.57831
M2 max thrust at cruise speed: 16.6682; M3 max thrust at cruise speed : 6.05204
Wing dimensions: 1.21626 x 0.334397
Passenger count: 32; minimum fairing floor area: 0.042366
Optimal battery capacity for M3 (for cruise only): 21.339; Max # of laps for M3: 6.04202 
*/
/*
    ofstream rec3("info_M3.txt");
    for (int i = 1; i <= runtime; i++) {
        x.randomize(part_num);
        if (max.score(mission) < x.score(mission)) {
            max.set(descent(x, 0.5, 0.5, 10000, mission, 1));
            for (int i = 0; i < x.status().size(); i++) {
                rec3 << max.status()[i] << ", ";
            }
            rec3 << max.score(mission) << endl;
        }

        for (int i = 0; i < x.status().size(); i++) {
            rec3 << x.status()[i] << ", ";
        }
        rec3 << x.score(mission) << endl;

        if (i % 10000 == 0) {
            max.paste(1);
        }
    }
    max.paste(1);
    rec3.close();

    if (max.score_M3() > m3_max) {
        m3_max = max.score_M3();
    }

    x.set(c);
    max.set(c);
    mission = 2;

    ofstream rec2("info_M2.txt");
    for (int i = 1; i <= runtime; i++) {
        x.randomize(part_num);
        if (max.score(mission) < x.score(mission)) {
            max.set(descent(x, 0.5, 0.5, 10000, mission, 1));
            for (int i = 0; i < x.status().size(); i++) {
                rec2 << max.status()[i] << ", ";
            }
            rec2 << max.score(mission) << endl;
        }

        for (int i = 0; i < x.status().size(); i++) {
            rec2 << x.status()[i] << ", ";
        }
        rec2 << x.score(mission) << endl;

        if (i % 10000 == 0) {
            max.paste(1);
        }
    }
    max.paste(1);
    rec2.close();

    if (max.score_M2() > m2_max) {
        m2_max = max.score_M2();
    }

    x.set(c);
    max.set(c);
    mission = 1;

    ofstream rec1("info_M1.txt");
    for (int i = 1; i <= runtime; i++) {
        x.randomize(part_num);
        if (max.score(mission) < x.score(mission)) {
            max.set(descent(x, 0.5, 0.5, 10000, mission, 1));
            for (int i = 0; i < x.status().size(); i++) {
                rec1 << max.status()[i] << ", ";
            }
            rec1 << max.score(mission) << endl;
        }

        for (int i = 0; i < x.status().size(); i++) {
            rec1 << x.status()[i] << ", ";
        }
        rec1 << x.score(mission) << endl;

        if (i % 10000 == 0) {
            max.paste(1);
        }
    }
    max.paste(1);
    rec1.close();
    max.paste(10);
*/