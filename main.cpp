/*
 * aige::SingingRivuleProject
 * file:main.cpp
 * GNU AFFERO GENERAL PUBLIC LICENSE
 * Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
 * Everyone is permitted to copy and distribute verbatim copies
 * of this license document, but changing it is not allowed.
 */
#include <soundtouch/SoundTouch.h>
#include <sox.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>

#include "freq.hpp"
#include "note.hpp"

using namespace std;

vmic::notes notes;
int sampleTime = 0;
int sampleRate = 0;

bool running = true;
int cmd = 0;

int main(int argc,const char *argv[]){
    assert(sox_init() == SOX_SUCCESS);
    sox_signalinfo_t s;
    s.channels = 1;
    s.length = 0;
    s.mult = 0;
    s.precision = 16;
    s.rate = 0;

    sox_format_t * in;
    sox_format_t * out;
    printf("open device\n");
    in  = sox_open_read("default",&s,NULL,"pulseaudio");
    if(!in){
        printf("fail to open device:default\n");
        exit(1);
    }
    const char * defaultOut = "null";
    const char * outName = defaultOut;
    if(argc>=2)
        outName = argv[1];
    out = sox_open_write(outName,&in->signal,NULL,"pulseaudio",NULL,NULL);
    if(!out){
        out = sox_open_write("default",&in->signal,NULL,"pulseaudio",NULL,NULL);
    }
    if(!out){
        printf("fail to open out device:%s\n",outName);
        exit(1);
    }

    signal(SIGINT,[](int){
        running = false;
    });
    signal(SIGUSR1,[](int){
        cmd = 1;
    });
    signal(SIGUSR2,[](int){
        cmd = 2;
    });

    auto st = new soundtouch::SoundTouch;
    st->setPitch(1);
    sampleRate = in->signal.rate;
    st->setSampleRate(sampleRate);
    st->setChannels(1);
    int32_t sampleBuffer_i[512];
    float   sampleBuffer_f[512];
    vmic::freqProcessor f;
    sampleTime = 0;
    notes.load(sampleRate);
    int lastTargetPitch = 0;

    while(running){

        switch (cmd) {
            case 1:
                printf("\n");
                notes.load(sampleRate);
                cmd = 0;
                break;
            case 2:
                printf("\nstart\n");
                sampleTime = 0;
                cmd = 0;
                break;
            default:
                {

                    size_t len = sox_read(in,sampleBuffer_i,512);
                    for(int i=0;i<512;++i){
                        sampleBuffer_f[i] = ((float)sampleBuffer_i[i])/((float)0x7fffffff);
                    }

                    int targetPitch = notes.getNote(sampleTime);
                    //如果目标改变了，清空频谱
                    if(targetPitch!=lastTargetPitch){
                        f.clearBuffer();
                    }
                    lastTargetPitch = targetPitch;
                    //获取频谱
                    f.pushBuffer(sampleBuffer_f);

                    for(int i=0;i<16;++i){
                        //设置变调目标
                        if(targetPitch>0){
                            f.getPitchShift(targetPitch);
                            st->setPitch(f.getPitchShiftCache());
                        }else{
                            st->setPitch(1);
                        }
                        st->putSamples(sampleBuffer_f+i*32 , 32);
                        st->receiveSamples(sampleBuffer_f+i*32 , 32);
                    }

                    for(int i=0;i<512;++i){
                        sampleBuffer_i[i] = sampleBuffer_f[i]*0x7fffffff;
                    }

                    printf("time=%-8d pitchNum=%-5lu target=%-5d shift=%-8f index=%-4d\r",
                           sampleTime/sampleRate,
                           f.pitches.size(),
                           targetPitch,
                           f.shiftOcatve,
                           notes.notes_index);
                    fflush(stdout);

                    sampleTime += len;
                    sox_write(out,sampleBuffer_i,len);
                }
                break;
        }
    }

    delete st;
    sox_close(out);
    sox_close(in);
    sox_quit();
    return 0;
}
