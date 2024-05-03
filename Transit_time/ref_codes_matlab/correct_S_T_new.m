function    [ S_T_new_2, out_q_rescaled]     = correct_S_T_new(S_T_old,out_q)

s_t_old                     = [S_T_old(1); diff(S_T_old)];
out_q_pdf                   = [out_q(1); diff(out_q)];
should_leave                = nansum(out_q_pdf);
% thresh                      = 10^-3;
thresh                      = 0;

% Zero pass
loc_neg_0                       = find(s_t_old <= thresh); 
loc_pos_0                       = find(s_t_old > thresh);
leftover_0                      = sum(out_q_pdf(loc_neg_0));
increment_0                     = leftover_0; 
out_q_pdf_rescaled_0            = out_q_pdf;
out_q_pdf_rescaled_0(loc_neg_0) = 0;
out_q_pdf_rescaled_0(loc_pos_0) = out_q_pdf_rescaled_0(loc_pos_0) + increment_0 / length(loc_pos_0);
check_0                         = nansum(out_q_pdf_rescaled_0) - should_leave;

% First pass
s_t_new_temp                    = s_t_old - out_q_pdf_rescaled_0;
loc_neg_1                       = find(s_t_new_temp < thresh); 
loc_pos_1                       = find(s_t_new_temp > thresh);
leftover_1                      = sum(s_t_new_temp(loc_neg_1));
out_q_pdf_rescaled_1            = out_q_pdf_rescaled_0;
out_q_pdf_rescaled_1(loc_neg_1) = s_t_old(loc_neg_1);
% s_t_new_temp(loc_neg_1)         = 0;                   
out_q_pdf_rescaled_1(loc_pos_1) = out_q_pdf_rescaled_1(loc_pos_1) - leftover_1/ length(loc_pos_1);
check_1                         = nansum(out_q_pdf_rescaled_1) - should_leave;

% Second pass
s_t_new_temp                    = s_t_old - out_q_pdf_rescaled_1;
loc_neg_2                       = find(s_t_new_temp < thresh); 
loc_pos_2                       = find(s_t_new_temp > thresh);
leftover_2                      = sum(s_t_new_temp(loc_neg_2));
out_q_pdf_rescaled_2            = out_q_pdf_rescaled_1;
out_q_pdf_rescaled_2(loc_neg_2) = s_t_old(loc_neg_2);
% s_t_new_temp(loc_neg_2)         = 0;                   
out_q_pdf_rescaled_2(loc_pos_2) = out_q_pdf_rescaled_2(loc_pos_2) - leftover_2/ length(loc_pos_2);
check_2                         = nansum(out_q_pdf_rescaled_2) - should_leave;


% % % % % figure(1)
% % % % % clf(1)
% % % % % plot(out_q_pdf_rescaled_0); hold on
% % % % % plot(out_q_pdf_rescaled_1); hold on
% % % % % plot(out_q_pdf_rescaled_2); hold on
% % % % % 
% % % % % xlim([ 0 1000])


out_q_rescaled              = cumsum(out_q_pdf_rescaled_2);
s_t_new_2                   = s_t_old - out_q_pdf_rescaled_2;
S_T_new_2                   = cumsum(s_t_new_2);

% has_to_leave                = out_q(end);
% left                        = -( nansum(s_t_new_2)-nansum(s_t_old) );
% remains                     = has_to_leave - left;
% diff_Q                      = nansum(out_q_pdf) - nansum(out_q_pdf_rescaled_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





