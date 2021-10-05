clear res


for subj = 1:length(d.subjects)
    res(subj).ee_rt = [];
    res(subj).he_rt = [];
    res(subj).eh_rt = [];
    res(subj).hh_rt = [];
    res(subj).ee_cor = [];
    res(subj).he_cor = [];
    res(subj).eh_cor = [];
    res(subj).hh_cor = [];
    
    trial = 0;
    while trial < length(d.subjects(subj).exp.rt)
        
        trial = trial+1;
        
        thisSubj = d.subjects(subj);
        
        if thisSubj.exp.stim_array{1,trial}.coh_difficulty == 1 && thisSubj.exp.stim_array{1,trial}.match_difficulty == 1
            res(subj).ee_rt = [res(subj).ee_rt,thisSubj.exp.rt(trial)];
            res(subj).ee_cor = [res(subj).ee_cor,thisSubj.exp.correct(trial)];
        elseif thisSubj.exp.stim_array{1,trial}.coh_difficulty == 1 && thisSubj.exp.stim_array{1,trial}.match_difficulty == 2
            res(subj).eh_rt = [res(subj).eh_rt,thisSubj.exp.rt(trial)];
            res(subj).eh_cor = [res(subj).eh_cor,thisSubj.exp.correct(trial)];
        elseif thisSubj.exp.stim_array{1,trial}.coh_difficulty == 2 && thisSubj.exp.stim_array{1,trial}.match_difficulty == 1
            res(subj).he_rt = [res(subj).he_rt,thisSubj.exp.rt(trial)];
            res(subj).he_cor = [res(subj).he_cor,thisSubj.exp.correct(trial)];
        elseif thisSubj.exp.stim_array{1,trial}.coh_difficulty == 2 && thisSubj.exp.stim_array{1,trial}.match_difficulty == 2
            res(subj).hh_rt = [res(subj).hh_rt,thisSubj.exp.rt(trial)];
            res(subj).hh_cor = [res(subj).hh_cor,thisSubj.exp.correct(trial)];
        end
        
    end
    
    res(subj).ee_meanrt = mean(res(subj).ee_rt,'omitnan');
    res(subj).ee_pc = accthis(res(subj).ee_cor);
    res(subj).eh_meanrt = mean(res(subj).eh_rt,'omitnan');
    res(subj).eh_pc = accthis(res(subj).eh_cor);
    res(subj).he_meanrt = mean(res(subj).he_rt,'omitnan');
    res(subj).he_pc = accthis(res(subj).he_cor);
    res(subj).hh_meanrt = mean(res(subj).hh_rt,'omitnan');
    res(subj).hh_pc = accthis(res(subj).hh_cor);
    
    
end

all.ee = mean(cell2mat({res(:).ee_meanrt}));
all.eh = mean(cell2mat({res(:).eh_meanrt}));
all.he = mean(cell2mat({res(:).he_meanrt}));
all.hh = mean(cell2mat({res(:).hh_meanrt}));
all.ee_pc = mean(cell2mat({res(:).ee_pc}));
all.eh_pc = mean(cell2mat({res(:).eh_pc}));
all.he_pc = mean(cell2mat({res(:).he_pc}));
all.hh_pc = mean(cell2mat({res(:).hh_pc}));

function accuracy = accthis(data)
accuracy = (sum(data)/length(data))*100;
end