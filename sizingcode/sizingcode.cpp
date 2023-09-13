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
