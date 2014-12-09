function asrObservations = runModelForward(p,A,B,C,D,v0,Ow,T,model)

asrObservations = zeros(T,1);
asrObservations(1) = v0;

knowledgeStates = zeros(T,1);

% TODO: why are rows flipped here? are they actually flipped?
A = flipud(A);
C = flipud(C);
D = flipud(D);

if strcmp(model, 'KTAttn')
    attentionStates = Ow;
    initialKnowledgeStateDist = p .* B(:,v0) .* D(:,attentionStates(1));
elseif strcmp(model, 'KT')
    initialKnowledgeStateDist = p .* B(:,v0);
else
    error('model %s not found!\n', model);
end
initialKnowledgeStateDist = initialKnowledgeStateDist./(sum(initialKnowledgeStateDist));
knowledgeStates(1) = discreteSample(initialKnowledgeStateDist,1);

for t=2:T
    
    if strcmp(model, 'KTAttn')
        
        % generate state
        previous_state = knowledgeStates(t-1);
        previous_attention_state = attentionStates(t-1);
        state_transition_probabilities = A(previous_state,:);
        attention_state_transition_probabilities = C(previous_attention_state,:);
        
        total_transition_prob = state_transition_probabilities.*attention_state_transition_probabilities;
        total_transition_prob = total_transition_prob./(sum(total_transition_prob));
        
        current_state = discreteSample(total_transition_prob,1);
        knowledgeStates(t) = current_state;
        current_attention_state = attentionStates(t);
        
        % generate observation for state
        performance_prob_for_knowledge_state = B(current_state,:);
        performance_prob_for_attention_state = D(current_attention_state,:);
        total_observation_prob = performance_prob_for_knowledge_state.*performance_prob_for_attention_state;
        total_observation_prob = total_observation_prob./(sum(total_observation_prob));
        asrObservations(t) = discreteSample(total_observation_prob,1);
        
        
    elseif strcmp(model, 'KT')
        
        % generate state
        previous_state = knowledgeStates(t-1);
        state_transition_probabilities = A(previous_state,:);  % TODO is this correct? seems like the rows are flipped
        current_state = discreteSample(state_transition_probabilities,1);
        knowledgeStates(t) = current_state;
        
        % generate observation for state
        performance_prob_for_knowledge_state = B(current_state,:);
        asrObservations(t) = discreteSample(performance_prob_for_knowledge_state,1);
        
    else
        error('model %s not found!\n', model);
    end
    
    
end
end