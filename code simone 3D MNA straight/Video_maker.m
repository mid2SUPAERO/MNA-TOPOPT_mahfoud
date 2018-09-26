clear
imageNames_conv = dir(fullfile('*convergence_*png'));
imageNames_conv = {imageNames_conv.name}';
imageNames_best = dir(fullfile('*configuration_*png'));
imageNames_best = {imageNames_best.name}';
jj=0;
jj1=0;
jj2=0;
outputVideo = VideoWriter(fullfile('convergence04.mp4'));
outputVideo.FrameRate = 10;
open(outputVideo)
for ii = 1:length(imageNames_conv)
    img1 = imread(fullfile(imageNames_conv{ii}));
    for ss=1:length(imageNames_best)
        if  strfind(imageNames_best{ss},sprintf('%04d',ii))
            jj=jj+1;
        end
    end
    img2 = imread(fullfile(imageNames_best{jj}));
    img=[img1,img2];
    writeVideo(outputVideo,img)
end

close(outputVideo)
% 
% % outputVideo.FrameRate = shuttleVideo.FrameRate;
% imageNames_animation = dir(fullfile('ani**ee.jpg'));
% imageNames_animation = {imageNames_animation.name}';
% outputVideo2 = VideoWriter(fullfile('animation.avi'));
% outputVideo2.FrameRate = 6;
% open(outputVideo2)
% for ii = 1:length(imageNames_animation)
%     img1 = imread(fullfile(imageNames_animation{ii}));
%     writeVideo(outputVideo2,img1)
% end
% for ii = 1:length(imageNames_animation)
%     img1 = imread(fullfile(imageNames_animation{ii}));
%     writeVideo(outputVideo2,img1)
% end
% close(outputVideo2)