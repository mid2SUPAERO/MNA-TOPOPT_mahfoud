clear
% imageNames_conv = dir(fullfile('*convergence_*png'));
% imageNames_conv = {imageNames_conv.name}';
imageNames_best = dir(fullfile('angle/*configuration_*png'));
imageNames_best = {imageNames_best.name}';
jj=0;
jj1=0;
jj2=0;
outputVideo = VideoWriter(fullfile('convergence05.avi'));
outputVideo.FrameRate = 12;
open(outputVideo)
for ii = 1:length(imageNames_best)
%     img1 = imread(fullfile(imageNames_conv{ii}));
%     for ss=1:length(imageNames_best)
%         if  strfind(imageNames_best{ss},sprintf('%03d',ii))
%             jj=jj+1;
%         end
%     end
    img2 = imread(fullfile(['angle/',imageNames_best{ii}]));
    img=[img2];
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