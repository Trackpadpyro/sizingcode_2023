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
        //t1.paste(state);
        switch (state) {
        case 1: score = t1.score();
            break;
        case 2: score = t1.score_M2();
            break;
        case 3: score = t1.score_M3();
            break;
        }
        if (j % 10 == 0 && j != 0) {
            if (temp3 == score) {
                break;
            }
            temp3 = score;
            if (j % 1000 == 0) {
                cout << j << ", " << temp3 << endl;
            }
        }
    }
    if (!silent) {
        cout << "End: " << endl;
        t1.paste(state);
    }
    else {
        cout << score << endl;
    }
    return t1;
}
int main()
{
    srand(time(NULL));


    int runtime = 1000000;
    vector<float> c = { 1.13293004, 0.492362887, 4.93311644, 8.18183994, 1.31787837, 415.073120, 18.2981548, 0.0776725709, 0.00242149946 };

    Parameter_Set test(l_bound);
    Parameter_Set max(c);
    int mission = 3;
    float score, score_temp;
    int counter = 0;

    max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
    ofstream file3("M3_info.txt");
    file3 << max.score_M3() << endl;

    while (true) {
        counter++;
        test.randomize(part_num);
        switch (mission) {
        case 1:
            score_temp = test.score();
            score = max.score();
            break;
        case 2:
            score_temp = test.score_M2();
            score = max.score_M2();
            break;
        case 3:
            score_temp = test.score_M3();
            score = max.score_M3();
            break;
        }
        if (score_temp > score) {
            max.set(test.status());
            max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
            file3 << max.score_M3() << endl;
            //max.paste(mission);
        }
        //cout << counter << ": " << score_temp  << ", " << score << endl;
        //cout << score_temp << endl;
        file3 << score_temp << endl;
        if (counter == runtime) {
            break;
        }
    }
    max.paste(1);
    file3.close();

    if (m3_max <= max.score_M3()) {
        m3_max = max.score_M3();
    }

    counter = 0;
    max.set(c);
    test.set(l_bound);
    score = 0;
    score_temp = 0;
    mission = 2;

    max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
    ofstream file2("M2_info.txt");
    file2 << max.score_M2() << endl;

    while (true) {
        counter++;
        test.randomize(part_num);
        switch (mission) {
        case 1:
            score_temp = test.score();
            score = max.score();
            break;
        case 2:
            score_temp = test.score_M2();
            score = max.score_M2();
            break;
        case 3:
            score_temp = test.score_M3();
            score = max.score_M3();
            break;
        }
        if (score_temp > score) {
            max.set(test.status());
            max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
            file3 << max.score_M2() << endl;
            //max.paste(mission);
        }
        //cout << counter << ": " << score_temp  << ", " << score << endl;
        //cout << score_temp << endl;
        file2 << score_temp << endl;
        if (counter == runtime) {
            break;
        }
    }
    max.paste(1);
    file2.close();

    if (m2_max <= max.score_M2()) {
        m2_max = max.score_M2();
    }

    counter = 0;
    max.set(c);
    test.set(l_bound);
    score = 0;
    score_temp = 0;
    mission = 1;

    max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
    ofstream file1("M_total_info.txt");
    file1 << max.score() << endl;

    while (true) {
        counter++;
        test.randomize(part_num);
        switch (mission) {
        case 1:
            score_temp = test.score();
            score = max.score();
            break;
        case 2:
            score_temp = test.score_M2();
            score = max.score_M2();
            break;
        case 3:
            score_temp = test.score_M3();
            score = max.score_M3();
            break;
        }
        if (score_temp > score) {
            max.set(test.status());
            max.set(descent(max, 0.5, 0.5, 10000, mission, 1));
            file1 << max.score() << endl;
            //max.paste(mission);
        }
        //cout << counter << ": " << score_temp  << ", " << score << endl;
        //cout << score_temp << endl;
        file1 << score_temp << endl;
        if (counter == runtime) {
            break;
        }
    }
    max.paste(1);
    file1.close();
}
