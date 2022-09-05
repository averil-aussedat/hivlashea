fr=10
rm *.mp4
ffmpeg -framerate 5 -i fi%d.png fi.mp4
ffmpeg -framerate $fr -i fe%d.png fe.mp4
ffmpeg -framerate $fr -i E%d.png E.mp4
ffmpeg -framerate $fr -i rho%d.png rho.mp4
