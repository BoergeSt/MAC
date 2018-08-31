function gcaExpandable()
    set(gca, 'ButtonDownFcn', [...
        'set(copyobj(gca, uipanel(''Position'', [0 0 1 1])), ' ...
        '    ''Units'', ''normal'', ''OuterPosition'', [0 0 1 1], ' ...
        '    ''ButtonDownFcn'', ''delete(get(gca, ''''Parent''''))''); ']);
    child_handles = allchild(gca);
    Turn_HitTest_off(child_handles)
end

function Turn_HitTest_off(child_handles)
    for x=child_handles
        set(x,'HitTest', 'off');
        %Turn_HitTest_off(allchild(x));
    end
end
