% combine evidence in each view 2022.3
function [multi_predict]=CombineEvidence(b,u,testnum,viewnum,im,class)
for jt=1:testnum
    S_a=0;
    for i=2:viewnum
        if S_a==0
            % b^(i-1) b^i.
            bb = b{i-1,im}(jt,:)' * b{i,im}(jt,:);
            % b^(i-1) * u^i
            bu = b{i-1,im}(jt,:) * u{i,im}(jt);
            % b^i * u^(i-1)
            ub = b{i,im}(jt,:) * u{i-1,im}(jt);
            % calculate C
            C = sum(sum(bb))-sum(diag(bb));

            % calculate b^a
            b_a = (b{i-1,im}(jt,:) .* b{i,im}(jt,:) + bu + ub) / (1-C);
            % calculate u^a
            u_a = u{i-1,im}(jt) * u{i,im}(jt)/(1-C);

            % calculate new S
            S_a = class / u_a;
            % calculate new e_k
            e_a = b_a * S_a;
            alpha_a = e_a + 1;
        else
            % b^(i-1) b^i
            bb = b_a' * b{i,im}(jt,:);
            % b^(i-1) * u^i
            bu = b_a * u{i,im}(jt);
            % b^i * u^(i-1)
            ub = b{i,im}(jt,:) * u_a;
            % calculate C
            C = sum(sum(bb))-sum(diag(bb));

            % calculate b^a
            b_a = (b_a .* b{i,im}(jt,:) + bu + ub) / (1-C);
            % calculate u^a
            u_a = u_a * u{i,im}(jt)/(1-C);

            % calculate new S
            S_a = class / u_a;
            % calculate new e_k
            e_a = b_a * S_a;
            alpha_a = e_a + 1;
        end
    end
    % final e_a
    [~, multi_e] = max(e_a);
    multi_predict{im}(jt) = multi_e;
end