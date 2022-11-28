[~, parcells] = getParcellsByLansAllansAtlas;
for i=1:size(parcells.names,1)
    parcells_new.indicators(:,:,i)=imtranslate(parcells.indicators(:,:,i),[0 -13]);
end
parcells_new.names=parcells.names;
parcells_new.description=parcells.description;
parcells_new.CombinedParcells=zeros(R,C);
exclude=[23,24,25,26]; %exclude sueprior and inferior colliculi in the combined image
for k=1:length(parcells_new.names)
    E=parcells_new.indicators(:,:,k)*k;
    if ismember(k,exclude)
        continue;
    end
    parcells_new.CombinedParcells=parcells_new.CombinedParcells+E;
end

%make a mask outsdie the brain and over the superior saggital sinus
%% mask the outside of the brain and superior saggital sinus along hte midline
mask=parcells_template>0;
figure;imshow(parcells_template); sinusMask=drawrectangle;sinusMaskidx=createMask(sinusMask);mask(sinusMaskidx)=0; imshow(mask);