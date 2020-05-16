function [flag,nr] = nonredundant(iid,tid,w)
    % Flags for redundant dummy levels.
    % Returns all the non-redundant levels of dummies
    obs = numel(iid);
    iid(w==0)=[];
    tid(w==0)=[];
    
    for k = 1:max(numel(iid),numel(tid))
        kid = where_id_with_single_obs(iid);
        iid(kid)=[];
        tid(kid)=[];
        kid = where_id_with_single_obs(tid);
        if ~any(kid), break; end
        iid(kid)=[];
        tid(kid)=[];
    end
    flag = obs ~= numel(iid);
    nr.iid = unique(iid);
    nr.tid = unique(tid);
    
    
    function kid = where_id_with_single_obs(id)
        kid = ~ismember(id,ids_with_multiple_obs(id));
    
    function kid = ids_with_multiple_obs(id)
        [~,wh]=unique(id);
        id(wh)=[];
        kid=unique(id);
        
        
        
        

                
    