
function [time_PLUS, freq_PLUS_new, da_dataPLUS_new, new_SCALE]=newTFRplot(oldTIME, oldFREQ, oldDATA, strechy)

da_dataPLUS=interp2(oldDATA,3);
time_PLUS=linspace(oldTIME(1),oldTIME(end),size(da_dataPLUS,2));
freq_PLUS=linspace(oldFREQ(1),oldFREQ(end),size(da_dataPLUS,1));
freq_PLUS_new=(linspace(freq_PLUS(1),freq_PLUS(end),floor(length(freq_PLUS))).^strechy)./100^(strechy-1);

new_index=[];
for iter=1:length(freq_PLUS_new)
[~,dummy]=min(abs(freq_PLUS_new(iter)-freq_PLUS));
new_index(iter)= dummy;
end
da_dataPLUS_new=da_dataPLUS(new_index,:);
new_SCALE = round(downsample(freq_PLUS_new,ceil(length(freq_PLUS_new)/freq_PLUS_new(end))));

end
















