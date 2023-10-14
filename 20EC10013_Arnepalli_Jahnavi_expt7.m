generatorset = [6 7 5];
stages = floor(log2(max(generatorset)));
num_states = 2^stages;
num_gen_pol = length(generatorset); %constraint length

% defining the encoder trellis as a Mealy Machine
% defining next states
for i=1:2.^stages
    for j=1:2
        next_states(i,j) = floor((i-1)/2) + (j-1)*(2.^(stages-1));
    end
end
% defining transition outputs
for i=1:2.^stages
    for j=1:2
        temp=0;
        for k=1:length(generatorset)
            tempval=rem(sum(bitand( transpose(int2bit(generatorset(k),stages+1)) , [j-1 transpose(int2bit((i-1),stages))])),2);
            temp = temp*2+tempval;
        end
        transition_outputs(i,j) = temp;
    end
end
% end of trellis creation

num_noises = 20;
noise = linspace(0,0.5,num_noises);
noise1 = linspace(0,0.5,num_noises);
bit_error_rate_p = zeros(num_noises,num_noises);

for(p = 1:1:num_noises)
    %p1 = p                      %used for BSC, for BAC comment this and uncomment the following loop
    for(p1 = 1:1:num_noises)
        for(iter=1:100)
            biterrorsrate = [];
            blockn = 8;     %the 256 bit message is broken into 16 blocks of length 16 each
            N = 32;
            for(block_iter=1:blockn)
                % Generate a vector of random bits using the randi() function
                d = randi([0 1], 1, N);

                % encoding sequence using above trellis - trellis is initialised to state 0
                state = 0;
                codedData = [];
                for i=1:N
                    codedData = cat(2,codedData,transpose(int2bit(transition_outputs(state+1, (d(i))+1),length(generatorset))));
                    state = next_states(state+1,(d(i))+1);
                end
                % transmit is the output from the convolutional encoder
                % verify the ratio of number of elements in data to transmit
                % that is the effective data rate
                BSCparameter1 = noise(p);
                BSCparameter2 = noise1(p1);
                codedData_org = codedData;
                for data=1:1:num_gen_pol*N
                    BSC = [];
                    if codedData_org(data) == 0
                        BSC = [BSC binornd(1,BSCparameter1,1,1)];
                    else
                        BSC = [BSC binornd(1,BSCparameter2,1,1)];
                    end
                end
                codedData = bitxor(codedData_org,BSC);
                states = [0];
                error = zeros(num_states,N);
                for i = 1 : 1 : num_states
                    for j = 1 : 1 : N
                        error(i, j) = 3*length(codedData);
                    end
                end

                for (i=1 : num_gen_pol : length(codedData)-(num_gen_pol-1))
                    new_states = [];
                    states = unique(states);
                    id  = floor( (i-1)/num_gen_pol) + 1;

                    for j= 1 : 1 : length(states)
                        index  =  states(j) + 1;
                        output_states = transition_outputs(index, :);
                        op1 = bitget(output_states(1), num_gen_pol : -1:1, 'int8');
                        op2 = bitget(output_states(2), num_gen_pol : -1:1, 'int8');
                        chunkData = codedData(1,i:i+(num_gen_pol-1));
                        error1 = 0;
                        error2 = 0;
                        for bit = 1 : 1: num_gen_pol
                            if(op1(bit)==0 && chunkData(bit)==1)
                                error1 = error1 + BSCparameter1*1.00;
                            end
                            if(op1(bit)==1 && chunkData(bit)==0)
                                error1 = error1 + BSCparameter2*1.00;
                            end
                            if(op1(bit)==0 && chunkData(bit)==0)
                                error1 = error1 + 1 - BSCparameter1*1.00;
                            end
                            if(op1(bit)==1 && chunkData(bit)==1)
                                error1 = error1 + 1 - BSCparameter2*1.00;
                            end
                        end
                        for bit = 1 : 1: num_gen_pol
                            if(op2(bit)==0 && chunkData(bit)==1)
                                error2 = error2 + BSCparameter1*1.00;
                            end
                            if(op2(bit)==1 && chunkData(bit)==0)
                                error2 = error2 + BSCparameter2*1.00;
                            end
                            if(op2(bit)==0 && chunkData(bit)==0)
                                error2 = error2 + 1- BSCparameter1*1.00;
                            end
                            if(op2(bit)==1 && chunkData(bit)==1)
                                error2 = error2 + 1 - BSCparameter2*1.00;
                            end
                        end
                        %disp(error1);
                        %disp(error2);
                        state1 = next_states(index,1);
                        state2 = next_states(index,2);
                        if id ==1
                            error(state1 + 1, id)= max(error(state1 + 1, id),error1);
                            error(state2 + 1, id)= max(error(state2 + 1, id),error2);
                        end
                        if id > 1
                            error(state1 + 1, id)= max(error(state1 + 1, id),error1 + error(index, id-1));
                            state2 = next_states(index,2);
                            error(state2 + 1, id)= max(error(state2 + 1, id),error2+ error(index, id-1));
                        end
                        new_states= [new_states, state1];
                        new_states= [new_states, state2];
                    end
                    states = new_states;
                end
                %%back tracking
                back_ptr = zeros(num_states, num_states);
                for i = 1 : 1: num_states
                    for j=1 : 1: num_states
                        back_ptr(i,j) = 0;
                    end
                end
                for i = 1 : 1: num_states
                    for j = 1 : 1 : length(next_states(i,:))
                        back_ptr(next_states(i,j)+1, i) = 1;
                    end
                end
                [last_cur_state, last_index] = max(error( :, N),[], "all");
                decoded_array = [];
                message_input = [];
                for i = N : -1: 1
                    back_error = N*num_gen_pol + 1;
                    back_state = last_index;
                    for j = 1 : 1 : num_states
                        if back_ptr(last_index, j)==1
                            if(i>1)
                                if back_error < error(j,i-1)
                                    back_error = max(back_error,error(j,i-1));
                                    back_state = j;
                                end
                            else
                                back_state = 1;
                            end
                        end
                    end

                    ptr1=[];
                    mark=zeros(2);

                    for j = 1 : 1 : num_states

                        if(next_states(j,1) + 1== last_index)
                            ptr1=[ptr1,j];
                            mark(length(ptr1))=1;
                        end
                        if(next_states(j,2) + 1== last_index)
                            ptr1 =[ptr1,j];
                            mark(length(ptr1))=2;
                        end
                    end

                    if(i>1)
                        if(error(ptr1(1),i-1)>error(ptr1(2),i-1))
                            bits = bitget(transition_outputs(ptr1(1),mark(1)), num_gen_pol : -1:1, 'int8');
                            decoded_array = [bits,decoded_array];
                            message_input = [mark(1)-1,message_input];
                        else
                            bits = bitget(transition_outputs(ptr1(2),mark(2)), num_gen_pol : -1:1, 'int8');
                            decoded_array = [bits,decoded_array];
                            message_input = [mark(2)-1,message_input];
                        end
                    else
                        bits = bitget(transition_outputs(1,mark(1)), num_gen_pol : -1:1, 'int8');
                        decoded_array = [bits,decoded_array];
                        message_input = [mark(1)-1,message_input];
                    end
                    last_index = back_state;
                end
                %disp(decoded_array)
                %disp(message_input);
                %disp(length(decoded_array));
                %disp(length(codedData));
                %disp(sum(bitxor(decoded_array,codedData_org)))
                %disp(sum(bitxor(d,message_input)))
                BER = (sum(bitxor(d,message_input))/N);
                biterrorsrate = [biterrorsrate BER];
            end
            bit_rate_p(iter) = mean(biterrorsrate);
        end
        bit_error_rate_p(p,p1) = mean(bit_rate_p);
    end
end

%plot BSCp BER
% for(i=1:num_noises)
%     bit_error_rate_BSCp(i) = bit_error_rate_p(i,i);
% end
% plot(noise,bit_error_rate_BSCp);
% xlabel("BSC parameter");
% ylabel("BER");

% surf() is used for 3D plotting in asymmetric case
surf(noise,noise1,bit_error_rate_p,'FaceAlpha',0.5);        %plot surface
