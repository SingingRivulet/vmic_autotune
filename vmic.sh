case $1 in
    start)
        echo 'start'
        killall vmic -12
    ;;
    load)
        echo 'load notes'
        if [  -n "$2" ]
        then
            echo 'downloading...'
            wget $2 -O notes.txt
        fi
        killall vmic -10
    ;;
    *)
        echo "使用方法："
        echo "重新播放：$0 start"
        echo "加载：$0 load [url]"
    ;;
esac
