/*
 * aige::SingingRivuleProject
 * file:freq.hpp
 * GNU AFFERO GENERAL PUBLIC LICENSE
 * Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 * Everyone is permitted to copy and distribute verbatim copies
 * of this license document, but changing it is not allowed.
 */
#ifndef FREQ_HPP
#define FREQ_HPP
#include <sox.h>
#include <vector>
#include <list>
#include <set>
#include <math.h>
#include <stdio.h>
#include <memory.h>

#define MIN_FREQ 20
namespace vmic{

struct cmplx{
    float r, i;
};

cmplx cmplx_mul_add(const cmplx c, const cmplx a, const cmplx b) {
    const cmplx ret = {
        (a.r * b.r) + c.r - (a.i * b.i),
        (a.i * b.r) + (a.r * b.i) + c.i
    };
    return ret;
}

void fft_Stockham(const cmplx *input, cmplx *output, int n, int flag) {
    int half = n >> 1;
    cmplx *buffer = (cmplx *) calloc(sizeof(cmplx), n * 2);
    if (buffer == NULL)
        return;
    cmplx *tmp = buffer;
    cmplx *y = tmp + n;
    memcpy(y, input, sizeof(cmplx) * n);
    for (int r = half, l = 1; r >= 1; r >>= 1) {
        cmplx *tp = y;
        y = tmp;
        tmp = tp;
        float factor_w = -flag * M_PI / l;
        cmplx w = {cosf(factor_w), sinf(factor_w)};
        cmplx wj = {1, 0};
        for (int j = 0; j < l; j++) {
            int jrs = j * (r << 1);
            for (int k = jrs, m = jrs >> 1; k < jrs + r; k++) {
                const cmplx t = {(wj.r * tmp[k + r].r) - (wj.i * tmp[k + r].i),
                                 (wj.i * tmp[k + r].r) + (wj.r * tmp[k + r].i)};
                y[m].r = tmp[k].r + t.r;
                y[m].i = tmp[k].i + t.i;
                y[m + half].r = tmp[k].r - t.r;
                y[m + half].i = tmp[k].i - t.i;
                m++;
            }
            const float t = wj.r;
            wj.r = (t * w.r) - (wj.i * w.i);
            wj.i = (wj.i * w.r) + (t * w.i);
        }
        l <<= 1;
    }
    memcpy(output, y, sizeof(cmplx) * n);
    free(buffer);
}

void fft_radix3(const cmplx *input, cmplx *output, int n, int flag) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int radix = 3;
    int np = n / radix;
    cmplx *res = (cmplx *) malloc(sizeof(cmplx) * n);
    cmplx *f0 = res;
    cmplx *f1 = f0 + np;
    cmplx *f2 = f1 + np;
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < radix; j++) {
            res[i + j * np] = input[radix * i + j];
        }
    }
    fft_radix3(f0, f0, np, flag);
    fft_radix3(f1, f1, np, flag);
    fft_radix3(f2, f2, np, flag);
    float wexp0 = -2 * M_PI * flag / n;
    cmplx wt = {cosf(wexp0), sinf(wexp0)};
    cmplx w0 = {1, 0};
    for (int i = 0; i < np; i++) {
        const float w0r = w0.r;
        w0.r = (w0r * wt.r) - (w0.i * wt.i);
        w0.i = (w0.i * wt.r) + (w0r * wt.i);
    }
    cmplx w = {1, 0};
    for (int j = 0; j < radix; j++) {
        cmplx wj = w;
        for (int k = 0; k < np; k++) {
            output[k + j * np] = cmplx_mul_add(f0[k], cmplx_mul_add(f1[k], f2[k], wj), wj);
            const float wjr = wj.r;
            wj.r = (wjr * wt.r) - (wj.i * wt.i);
            wj.i = (wj.i * wt.r) + (wjr * wt.i);
        }
        const float wr = w.r;
        w.r = (wr * w0.r) - (w.i * w0.i);
        w.i = (w.i * w0.r) + (wr * w0.i);
    }
    free(res);
}

void fft_radix5(const cmplx *input, cmplx *output, int n, int flag) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int radix = 5;
    int np = n / radix;
    cmplx *res = (cmplx *) calloc(sizeof(cmplx), n);
    cmplx *f0 = res;
    cmplx *f1 = f0 + np;
    cmplx *f2 = f1 + np;
    cmplx *f3 = f2 + np;
    cmplx *f4 = f3 + np;
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < radix; j++) {
            res[i + j * np] = input[radix * i + j];
        }
    }
    fft_radix5(f0, f0, np, flag);
    fft_radix5(f1, f1, np, flag);
    fft_radix5(f2, f2, np, flag);
    fft_radix5(f3, f3, np, flag);
    fft_radix5(f4, f4, np, flag);
    float wexp0 = -2 * M_PI * flag / n;
    cmplx wt = {cosf(wexp0), sinf(wexp0)};
    cmplx w0 = {1, 0};
    for (int i = 0; i < np; i++) {
        const float w0r = w0.r;
        w0.r = (w0r * wt.r) - (w0.i * wt.i);
        w0.i = (w0.i * wt.r) + (w0r * wt.i);
    }
    cmplx w = {1, 0};
    for (int j = 0; j < radix; j++) {
        cmplx wj = w;
        for (int k = 0; k < np; k++) {
            output[k + j * np] = cmplx_mul_add(f0[k], cmplx_mul_add(f1[k], cmplx_mul_add(f2[k],
                                                                                         cmplx_mul_add(f3[k], f4[k],
                                                                                                       wj), wj), wj),
                                               wj);
            const float wjr = wj.r;
            wj.r = (wjr * wt.r) - (wj.i * wt.i);
            wj.i = (wj.i * wt.r) + (wjr * wt.i);
        }
        const float wr = w.r;
        w.r = (wr * w0.r) - (w.i * w0.i);
        w.i = (w.i * w0.r) + (wr * w0.i);
    }
    free(res);
}

void fft_radix6(const cmplx *input, cmplx *output, int n, int flag) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int radix = 6;
    int np = n / radix;
    cmplx *res = (cmplx *) calloc(sizeof(cmplx), n);
    cmplx *f0 = res;
    cmplx *f1 = f0 + np;
    cmplx *f2 = f1 + np;
    cmplx *f3 = f2 + np;
    cmplx *f4 = f3 + np;
    cmplx *f5 = f4 + np;
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < radix; j++) {
            res[i + j * np] = input[radix * i + j];
        }
    }
    fft_radix6(f0, f0, np, flag);
    fft_radix6(f1, f1, np, flag);
    fft_radix6(f2, f2, np, flag);
    fft_radix6(f3, f3, np, flag);
    fft_radix6(f4, f4, np, flag);
    fft_radix6(f5, f5, np, flag);
    float wexp0 = -2 * M_PI * flag / n;
    cmplx wt = {cosf(wexp0), sinf(wexp0)};
    cmplx w0 = {1, 0};
    for (int i = 0; i < np; i++) {
        const float w0r = w0.r;
        w0.r = (w0r * wt.r) - (w0.i * wt.i);
        w0.i = (w0.i * wt.r) + (w0r * wt.i);
    }
    cmplx w = {1, 0};
    for (int j = 0; j < radix; j++) {
        cmplx wj = w;
        for (int k = 0; k < np; k++) {
            output[k + j * np] = cmplx_mul_add(f0[k], cmplx_mul_add(f1[k], cmplx_mul_add(f2[k],
                                                                                         cmplx_mul_add(f3[k],
                                                                                                       cmplx_mul_add(
                                                                                                           f4[k],
                                                                                                           f5[k],
                                                                                                           wj), wj),
                                                                                         wj), wj), wj);
            const float wjr = wj.r;
            wj.r = (wjr * wt.r) - (wj.i * wt.i);
            wj.i = (wj.i * wt.r) + (wjr * wt.i);
        }
        const float wr = w.r;
        w.r = (wr * w0.r) - (w.i * w0.i);
        w.i = (w.i * w0.r) + (wr * w0.i);
    }
    free(res);
}

void fft_radix7(const cmplx *input, cmplx *output, int n, int flag) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int radix = 7;
    int np = n / radix;
    cmplx *res = (cmplx *) calloc(sizeof(cmplx), n);
    cmplx *f0 = res;
    cmplx *f1 = f0 + np;
    cmplx *f2 = f1 + np;
    cmplx *f3 = f2 + np;
    cmplx *f4 = f3 + np;
    cmplx *f5 = f4 + np;
    cmplx *f6 = f5 + np;
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < radix; j++) {
            res[i + j * np] = input[radix * i + j];
        }
    }
    fft_radix7(f0, f0, np, flag);
    fft_radix7(f1, f1, np, flag);
    fft_radix7(f2, f2, np, flag);
    fft_radix7(f3, f3, np, flag);
    fft_radix7(f4, f4, np, flag);
    fft_radix7(f5, f5, np, flag);
    fft_radix7(f6, f6, np, flag);
    float wexp0 = -2 * M_PI * flag / n;
    cmplx wt = {cosf(wexp0), sinf(wexp0)};
    cmplx w0 = {1, 0};
    for (int i = 0; i < np; i++) {
        const float w0r = w0.r;
        w0.r = (w0r * wt.r) - (w0.i * wt.i);
        w0.i = (w0.i * wt.r) + (w0r * wt.i);
    }
    cmplx w = {1, 0};
    for (int j = 0; j < radix; j++) {
        cmplx wj = w;
        for (int k = 0; k < np; k++) {
            output[k + j * np] = cmplx_mul_add(f0[k], cmplx_mul_add(f1[k], cmplx_mul_add(f2[k],
                                                                                         cmplx_mul_add(f3[k],
                                                                                                       cmplx_mul_add(
                                                                                                           f4[k],
                                                                                                           cmplx_mul_add(
                                                                                                               f5[k],
                                                                                                               f6[k],
                                                                                                               wj),
                                                                                                           wj), wj),
                                                                                         wj), wj), wj);
            const float wjr = wj.r;
            wj.r = (wjr * wt.r) - (wj.i * wt.i);
            wj.i = (wj.i * wt.r) + (wjr * wt.i);
        }
        const float wr = w.r;
        w.r = (wr * w0.r) - (w.i * w0.i);
        w.i = (w.i * w0.r) + (wr * w0.i);
    }
    free(res);
}

void fft_Bluestein(const cmplx *input, cmplx *output, int n, int flag) {
    int m = 1 << ((unsigned int) (ilogbf((float) (2 * n - 1))));
    if (m < 2 * n - 1) {
        m <<= 1;
    }
    cmplx *y = (cmplx *) calloc(sizeof(cmplx), 3 * m);
    if (y == NULL)
        return;
    cmplx *w = y + m;
    cmplx *ww = w + m;
    w[0].r = 1;
    if (flag == -1) {
        y[0].r = input[0].r;
        y[0].i = -input[0].i;
        for (int i = 1; i < n; i++) {
            const float wexp = M_PI * i * i / n;
            w[i].r = cosf(wexp);
            w[i].i = sinf(wexp);
            w[m - i] = w[i];
            y[i].r = (input[i].r * w[i].r) - (input[i].i * w[i].i);
            y[i].i = (-input[i].i * w[i].r) - (input[i].r * w[i].i);
        }
    } else {
        y[0].r = input[0].r;
        y[0].i = input[0].i;
        for (int i = 1; i < n; i++) {
            const float wexp = M_PI * i * i / n;
            w[i].r = cosf(wexp);
            w[i].i = sinf(wexp);
            w[m - i] = w[i];
            y[i].r = (input[i].r * w[i].r) + (input[i].i * w[i].i);
            y[i].i = (input[i].i * w[i].r) - (input[i].r * w[i].i);
        }
    }
    fft_Stockham(y, y, m, 1);
    fft_Stockham(w, ww, m, 1);
    for (int i = 0; i < m; i++) {
        const float r = y[i].r;
        y[i].r = (r * ww[i].r) - (y[i].i * ww[i].i);
        y[i].i = (y[i].i * ww[i].r) + (r * ww[i].i);
    }
    fft_Stockham(y, y, m, -1);
    if (flag == -1) {
        for (int i = 0; i < n; i++) {
            output[i].r = ((y[i].r * w[i].r) + (y[i].i * w[i].i)) / m;
            output[i].i = -((y[i].i * w[i].r) - (y[i].r * w[i].i)) / m;
        }
    } else {
        for (int i = 0; i < n; i++) {
            output[i].r = ((y[i].r * w[i].r) + (y[i].i * w[i].i)) / m;
            output[i].i = ((y[i].i * w[i].r) - (y[i].r * w[i].i)) / m;
        }
    }
    free(y);
}


int base(int n) {
    int t = n & (n - 1);
    if (t == 0) {
        return 2;
    }
    for (int i = 3; i <= 7; i++) {
        int n2 = n;
        while (n2 % i == 0) {
            n2 /= i;
        }
        if (n2 == 1) {
            return i;
        }
    }
    return n;
}

void FFT(const cmplx *input, cmplx *output, int n) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int p = base(n);
    switch (p) {
        case 2:
            fft_Stockham(input, output, n, 1);
            break;
        case 3:
            fft_radix3(input, output, n, 1);
            break;
        case 5:
            fft_radix5(input, output, n, 1);
            break;
        case 6:
            fft_radix6(input, output, n, 1);
            break;
        case 7:
            fft_radix7(input, output, n, 1);
            break;
        default:
            fft_Bluestein(input, output, n, 1);
            break;
    }
}

void IFFT(const cmplx *input, cmplx *output, int n) {
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    int p = base(n);
    switch (p) {
        case 2:
            fft_Stockham(input, output, n, -1);
            break;
        case 3:
            fft_radix3(input, output, n, -1);
            break;
        case 5:
            fft_radix5(input, output, n, -1);
            break;
        case 6:
            fft_radix6(input, output, n, -1);
            break;
        case 7:
            fft_radix7(input, output, n, -1);
            break;
        default: {
                fft_Bluestein(input, output, n, -1);
                break;
            }
    }
    for (int i = 0; i < n; i++) {
        output[i].r = output[i].r / n;
        output[i].i = output[i].i / n;
    }
}

class freqProcessor{
    public:
        cmplx wave[8192];
        cmplx freq[8192];
        float freq_r[8192];
        float window[8192];
        std::list<float*> buffer;
        std::set<int> pitches;

        int pitch;
        float pitchShift;

        freqProcessor(){
            for(int i=0;i<16;++i){
                buffer.push_back(new float[512]);
            }
            pitchShift = 1;
            for(int i=0;i<8192;++i) {
                window[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (8192.0 - 1.0)));
            }
        }
        ~freqProcessor(){
            for(auto it:buffer){
                delete [] it;
            }
        }
        void pushBuffer(const float * in){
            auto t = buffer.front();
            buffer.pop_front();
            for(int i=0;i<512;++i){
                t[i] = in[i];
            }
            buffer.push_back(t);
            this->process();
        }
        void process(){
            int index = 0;
            for(auto & it:buffer){
                for(int i=0;i<512;++i){
                    wave[index].r = it[i]*window[index];
                    wave[index].i = 0;
                    ++index;
                }
            }
            FFT(wave,freq,8192);
            buildRealFreq();
            getPitch();
        }
        void buildRealFreq(){
            for(int i=0;i<8192;++i){
                freq_r[i] = freq[i].i*freq[i].i + freq[i].r*freq[i].r;
            }
        }
        void getPitch(){
            int i,value,maxValue=0,maxPosi=0;
            pitches.clear();
            for(i=0;i<600;++i){
                value = freq_r[i];
                if(i > MIN_FREQ && i < 600 && value>64*64 &&
                   value > freq_r[i - 1] &&
                   value >= freq_r[i + 1] &&
                   value > freq_r[i - 2] &&
                   value >= freq_r[i + 2] &&
                   value > freq_r[i - 3] &&
                   value >= freq_r[i + 3]) {
                    pitches.insert(i);
                    if(value > maxValue){
                        maxValue = value;
                        maxPosi = i;
                    }
                }
            }
            pitch = maxPosi;
            //printf("\n%d\n",pitches.size());
        }
#define log2 0.30102999566398
        int getNearPitch(int f1){
            if(!pitches.empty()){
                auto it = pitches.upper_bound(f1);//获取大于的第一个元素
                if(it==pitches.end()){//指向结尾
                    it--;
                    return *it;//返回最后一个元素
                }
                if(it==pitches.begin()){//指向开头
                    return *it;
                }
                int h = *it;
                it--;
                int l = *it;
                float h_d = fabs(log((float)h / (float)f1) / log2);
                float l_d = fabs(log((float)l / (float)f1) / log2);
                if(h_d<l_d){
                    return h;
                }else{
                    return l;
                }
            }
            return -1;
        }
        float getPitchShift(int f1){
            if(pitch>MIN_FREQ){
                /*
                int f[3];
                //上下八度查找
                f[0]=f1;
                f[1]=f1>>1;
                f[2]=f1<<1;
                float min_d = fabs(log((float)pitch / (float)f1) / log2);
                int min_i = 0;
                for(int i=1;i<3;++i){
                    float d = fabs(log((float)pitch / (float)f[i]) / log2);
                    if(d<min_d){
                        min_d = d;
                        min_i = i;
                    }
                }
                float targf = f[min_i];
                pitchShift = (((float)targf)/((float)pitch));
                return pitchShift;
                */
                int p = getNearPitch(f1);
                if(p>0){
                    pitchShift = (((float)f1)/((float)p));
                }
            }
            return this->pitchShift;
        }
        void clearBuffer(){
            for(auto & it:buffer){
                for(int i=0;i<512;++i){
                    it[i] = 0;
                }
            }
        }
};

}
#endif // FREQ_HPP
