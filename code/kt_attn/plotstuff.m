
figure; 
    subplot(1,2,1); hold on;
        plot(repmat(p_learn_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_learn_KT,'ro', 'LineWidth', 3, 'MarkerSize', 5)
        plot(all_p_learn_KTAttn,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Learn Rate (l)')
        ylim([-.1 1]);
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        ylim([-.1 1]);
        
    subplot(1,2,2); hold on;
        plot(repmat(p_forget_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_forget_KT,'ro', 'LineWidth', 3, 'MarkerSize', 5)
        plot(all_p_forget_KTAttn,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Forget Rate (f)')
        ylim([-.1 1]);
        xlabel('Subject ID')
        legend('Seed', 'KT', 'KT-Attn')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)

    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);

    
    
    
    %%
    
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_guess_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_guess_KT,'ro', 'LineWidth', 3, 'MarkerSize', 5);
        plot(all_p_guess_KTAttn,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Guess Rate (g)')
        ylim([-.1 1.1]);
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        
    subplot(1,2,2); hold on;
        plot(repmat(p_slip_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_slip_KT,'ro', 'LineWidth', 3, 'MarkerSize', 5);
        plot(all_p_slip_KTAttn,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Slip Rate (s)')
        ylim([-.1 1.1]);
        legend('Seed', 'KT', 'KT-Attn')
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)

    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);

    
    %%
    
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_learn_a_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_learn_a,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Attention Learn Rate (l_a)')
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        ylim([-.1 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_dontlearn_a_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_dontlearn_a,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Attention No-Learn Rate (f_a)')
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        ylim([-.1 1]);
        legend('Seed', 'KT-Attn')
    
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);
    
figure;
    subplot(1,2,1); hold on;
        plot(repmat(p_guess_a_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_guess_a,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Attention Guess Rate (g_a)')
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        ylim([-.1 1]);
    subplot(1,2,2); hold on;
        plot(repmat(p_slip_a_0, N, 1), 'k--', 'LineWidth', 3);
        plot(all_p_slip_a,'bo', 'LineWidth', 3, 'MarkerSize', 5)
        title('Attention Slip Rate (s_a)')
        xlabel('Subject ID')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',14)
        ylim([-.1 1]);
        legend('Seed', 'KT-Attn')
    
    
    figureHandle = gcf;
    set(findall(figureHandle,'type','text'),'fontSize',14);
    
    
    %%

total_num_sequences = 0;
    
load('subjects.mat');
N = length(subjectids);
for i=1:N
    
    seqs_filename = char(strcat('subjects/', subjectid, '_sequences.mat'));
    load(seqs_filename);
    sN = length(sequences.accept);
    
    total_num_sequences = total_num_sequences + sN;
end

total_num_sequences