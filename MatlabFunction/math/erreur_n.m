function [ erreur,norm_e ] = erreur_n( func_target, func2 )

%% do normalise error
%% [ erreur,norm_e ] = erreur_n( func_target, func2 )
erreur=sqrt(abs(func_target-func2).^2./abs(func_target).^2)*100;
norm_e=norm(func_target-func2)./norm(func_target);


end

