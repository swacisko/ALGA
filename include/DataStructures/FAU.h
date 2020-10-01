//
// Created by sylwester on 1/14/19.
//

#ifndef GENOMEALIGNMENT_FAU_H
#define GENOMEALIGNMENT_FAU_H


class FAU {
public:
    int *p, *w;

    FAU(int n) : p(new int[n]), w(new int[n]) {
        for (int x = 0; x < n; x++) {
            p[x] = -1;
            w[x] = 1;
        }
    }

    ~FAU() {
        delete[] p;
        delete[] w;
    }

    int Find(int x) {
        return (p[x] < 0) ? x : p[x] = Find(p[x]);
    }

    int Weight(int x) {
        return w[Find(x)];
    }

    void Union(int x, int y) {
        if ((x = Find(x)) == (y = Find(y))) return;
        if (w[x] > w[y]) {
            p[y] = x;
            w[x] += w[y];
        } else {
            p[x] = y;
            w[y] += w[x];
        }


    }

};


#endif //GENOMEALIGNMENT_FAU_H
