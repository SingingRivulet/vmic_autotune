#ifndef NOTE_HPP
#define NOTE_HPP
#include <vector>
#include <tuple>
#include <string>

namespace vmic{

class notes{
    public:
        std::vector<std::tuple<int,int,int> > note;
        int notes_index;

        notes(){
            notes_index = 0;
        }

        void load(int sampleRate,const std::string & path = "./notes.txt"){
            printf("load notes\n");
            FILE * fp = fopen(path.c_str(),"r");
            char buf[4096];
            float time_start,time_dur;
            int pitch;
            if(fp){
                note.clear();
                fgets(buf,sizeof(buf),fp);
                int notes_num = atoi(buf);
                if(notes_num>0){
                    printf("notes num:%d\n",notes_num);
                    for(int i=0;(i<notes_num && !feof(fp));++i){
                        bzero(buf,sizeof(buf));
                        fgets(buf,sizeof(buf),fp);
                        if(sscanf(buf,"%f %f %d",&time_start,&time_dur,&pitch)>=3){
                            int start = time_start*sampleRate;
                            int dur   = time_dur*sampleRate;
                            int targetPitch = (440 * 8192 * pow(2 , (pitch-57)/12 ))/sampleRate;
                            note.push_back(std::make_tuple(start,dur,targetPitch));
                        }
                    }
                }
                fclose(fp);
            }
        }

        int inArea(const std::tuple<int,int,int> & n,int t){
            if(t>=std::get<0>(n)){
                if(t<(std::get<0>(n)+std::get<1>(n))){
                    return 0;
                }else{
                    return 1;
                }
            }else{
                return -1;
            }
        }

        int getNoteIndex(int time){
            int status;
            if(note.empty()){
                return -1;
            }else{

                if(notes_index>=(int)note.size()){
                    notes_index = note.size()-1;
                }else if(notes_index<0){
                    notes_index = 0;
                }

                status = inArea(note[notes_index],time);//检查当前位置
                if(status==0){
                    return notes_index;
                }else if(status==1){
                    for( ; notes_index<(int)note.size() ; ++notes_index){
                        if(inArea(note[notes_index],time)==0){
                            return notes_index;
                        }
                    }
                }else if(status==-1){
                    for( ; notes_index>=0 ; --notes_index){
                        if(inArea(note[notes_index],time)==0){
                            return notes_index;
                        }
                    }
                }
            }
            return -1;
        }

        int getNote(int time){
            int id = getNoteIndex(time);
            try{
                if(id>=0){
                    return std::get<2>(note.at(id));
                }
            }catch(...){}
            return 0;
        }
};

}
#endif // NOTE_HPP
