function[h1,h2,handl1,handl2,Color_data]= plotPopAverage(h1,h2,Color_data,eventName, NormType,BrainReg,fsimaging,preEventWin, postEventWin,ylab,ylimit,color,figurecolor)
%plot mean +/-sem plots
% mean and sem across animals
numAnimals=size(Color_data.(eventName).(NormType).(BrainReg),2);
Color_data.(eventName).meanPop.(NormType).(BrainReg)=nanmean(Color_data.(eventName).(NormType).(BrainReg),2);
Color_data.(eventName).stdPop.(NormType).(BrainReg)=nanstd(Color_data.(eventName).(NormType).(BrainReg),0,2);
Color_data.(eventName).semPop.(NormType).(BrainReg)=Color_data.(eventName).stdPop.(NormType).(BrainReg)./sqrt(numAnimals);
upperSEM_animals=Color_data.(eventName).meanPop.(NormType).(BrainReg)+Color_data.(eventName).semPop.(NormType).(BrainReg);
lowerSEM_animals=Color_data.(eventName).meanPop.(NormType).(BrainReg)-Color_data.(eventName).semPop.(NormType).(BrainReg);
%mean and sem across trials of example mouse
numTrials=size(Color_data.(eventName).(NormType).(strcat(BrainReg,'_ExampleMouse')),2);
Color_data.(eventName).meanPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))=nanmean(Color_data.(eventName).(NormType).(strcat(BrainReg,'_ExampleMouse')),2);
Color_data.(eventName).stdPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))=nanstd(Color_data.(eventName).(NormType).(strcat(BrainReg,'_ExampleMouse')),0,2);
Color_data.(eventName).semPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))=Color_data.(eventName).stdPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))./sqrt(numTrials);
upperSEM_extrials=Color_data.(eventName).meanPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))+Color_data.(eventName).semPop.(NormType).(strcat(BrainReg,'_ExampleMouse'));
lowerSEM_extrials=Color_data.(eventName).meanPop.(NormType).(strcat(BrainReg,'_ExampleMouse'))-Color_data.(eventName).semPop.(NormType).(strcat(BrainReg,'_ExampleMouse'));

%plot figures;
%plot average +/- sem across animals
tstamps=((-preEventWin+1/fsimaging):1/fsimaging:(postEventWin))';
set(0, 'CurrentFigure', h1);handl1=plot(tstamps,Color_data.(eventName).meanPop.(NormType).(BrainReg),'Color',figurecolor,'LineWidth',2);box off;hold on;
xfill=[tstamps; flipud(tstamps)];
yfill=[lowerSEM_animals;flipud(upperSEM_animals)];
fill(xfill,yfill,'k','EdgeColor','None','facealpha',.5) ;
title(strcat(color,'-',eventName,'-Animals')); ylim(ylimit); ylabel(ylab);xlabel('Time(s)');

%plot average +/- sem across trials of example mouse
set(0, 'CurrentFigure', h2);handl2=plot(tstamps,Color_data.(eventName).meanPop.(NormType).(strcat(BrainReg,'_ExampleMouse')),'Color',figurecolor,'LineWidth',2);box off;hold on;
xfill=[tstamps; flipud(tstamps)];
yfill=[lowerSEM_extrials;flipud(upperSEM_extrials)];
fill(xfill,yfill,'k','EdgeColor','None','facealpha',.5) ;
title(strcat(color,'-',eventName,'-ExampleAnimalTrials')); ylim(ylimit); ylabel(ylab);xlabel('Time(s)');
end

