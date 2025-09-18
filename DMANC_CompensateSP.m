%% Distributed MANC depending on gradient tranmission
% 1 reference, K secondary source, K error sensors


classdef DMANC_CompensateSP
    properties
        Wc     % global control filter (K x (wlen+clen-1))
        Nabla  % local gradient (K x wlen)
        wlen   % length of control filter
        SecP   % secondary path estimates (K x K x slen)
        slen  % length of secondary path
        Xd    % Disturbance
        yc    % control signal
        Numnode  % number of nodes, K
        C     % compensate filter (M x M x clen)
        clen  % length of compensate filter
    end


    methods
        % initial
        function obj = DMANC_CompensateSP(wLen,SecondaryPath,sLen,node_num,N,Dis,cLen)
            obj.wlen    = wLen;
            obj.Wc      = zeros(node_num,(wLen));
            obj.Nabla   = zeros(node_num,(wLen+cLen-1));
            obj.SecP    = SecondaryPath;
            obj.slen    = sLen;
            obj.Numnode = node_num;
            obj.yc      = zeros(node_num,N);
%             obj.Xd      = zeros(node_num,N);
            obj.Xd = Dis;
%             for i = 1:node_num
%                 obj.Xd(i,:) = filter(PrimaryPath(i,:),1,xin);
%             end
            obj.C       = zeros(node_num,node_num,cLen);
            obj.clen    = cLen;
        end

       % obtain compensation filter
        function [err,obj] = CompensateSP(obj,muc)
            T   = 200000;           % duration
            wgn = randn(1,T);       % generate white noise
            err = zeros(obj.Numnode,obj.Numnode,T); % error signal

            for m = 1: obj.Numnode
                for k = 1: obj.Numnode
                    if m == k
                        continue;
                    else
                        wgnd = filter(reshape(obj.SecP(m,k,:),[1,obj.slen]),1,wgn);                     % wgn pass cross secondary path
                        FwgnLMS = dsp.FilteredXLMSFilter(obj.clen,"StepSize",muc,...
                        "SecondaryPathCoefficients",reshape(obj.SecP(m,m,:),[1,obj.slen]),...
                        "SecondaryPathEstimate",reshape(obj.SecP(m,m,:),[1,obj.slen]));
                        [~,e] = FwgnLMS(wgn,wgnd);
                        err(m,k,:) = e;
                        CF = -flip(FwgnLMS.Coefficients);
                        obj.C(m,k,:) = reshape(CF,[1,1,obj.clen]);
                    end
                end
            end

        end

        % normal no delay; transmit gradient; 144
        function [e,obj] = DMANC_gradient(obj,xin,muw)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            for i =1:N
                xc = [xin(i) xc(1:(end-1))];   % update reference vector
                

                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));
                e(4,i) = obj.Xd(4,i)-y4;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = muw*xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = muw*xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = muw*xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = muw*xf(4,:)*e(4,i);


                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);

            end


        end

        % normal no delay; transmit gradient; 166
        function [e,obj] = DMANC_gradient_166(obj,xin,muw)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            for i =1:N
                xc = [xin(i) xc(1:(end-1))];   % update reference vector
                

                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = muw*xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = muw*xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = muw*xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = muw*xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = muw*xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = muw*xf(6,:)*e(6,i);


                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

            end


        end

        %  delay; transmit gradient
        function [e,Delay_time,MUW,obj] = DMANC_gradient_vss(obj,xin,muw,delay_time,maxdelay)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(1,N);
            Delay_time = zeros(1,N + maxdelay);
%             Delay = delay_time;


            for i =1:N

                Delay = round(rand(1)*16000);
                Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
                muw = 1e-5*exp(-1*(Delay_time(i)/16000));
                MUW(i) = muw;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay,:) = Nabla1_store(i+Delay,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay,:) = Nabla2_store(i+Delay,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay,:) = Nabla3_store(i+Delay,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay,:) = Nabla4_store(i+Delay,:) + obj.Nabla(4,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end);


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));
                e(4,i) = obj.Xd(4,i)-y4;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = muw*xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = muw*xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = muw*xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = muw*xf(4,:)*e(4,i);


            end
        end


        %  delay; transmit gradient; 166 ; case 4; fluctuating network
        function [e,Delay_time,Delay_store,MUW,obj] = DMANC_gradient_vss_case4(obj,xin,muw0,maxdelay)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(1,N);
            Delay_time = zeros(1,N + maxdelay);
            Delay_store = zeros(1,N);
%             Delay = delay_time;


            for i =1:N

%                 Delay = round(rand(1)*16000);
                Delay = round((sin(2*pi*0.1*(i/16000)-(pi/2))+1)*8000);
                Delay_store(i) = Delay; 
                Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
                if Delay ~= 0 && Delay_time(i) == 0 && i > 1
                    Delay_time(i) = Delay_time(i-1);
                end
                muw = muw0*exp(-2*(Delay_time(i)/16000));
                MUW(i) = muw;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay,:) = Nabla1_store(i+Delay,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay,:) = Nabla2_store(i+Delay,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay,:) = Nabla3_store(i+Delay,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay,:) = Nabla4_store(i+Delay,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay,:) = Nabla5_store(i+Delay,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay,:) = Nabla6_store(i+Delay,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end

        %  delay; transmit gradient; 166 ; case 4; no vss
        function [e,Delay_time,MUW,obj] = DMANC_gradient_novss_case4(obj,xin,muw0,maxdelay,delaystore)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(1,N);
            Delay_time = zeros(1,N + maxdelay);
%             Delay = delay_time;


            for i =1:N

%                 Delay = round(rand(1)*16000);
                Delay = delaystore(i);
                Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
                muw = muw0;
                MUW(i) = muw;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay,:) = Nabla1_store(i+Delay,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay,:) = Nabla2_store(i+Delay,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay,:) = Nabla3_store(i+Delay,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay,:) = Nabla4_store(i+Delay,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay,:) = Nabla5_store(i+Delay,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay,:) = Nabla6_store(i+Delay,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end
    
        %  delay; transmit gradient; 166 ; case 3; sudden change delay
        function [e,Delay_time,MUW,obj] = DMANC_gradient_vss_case3(obj,xin,muw0)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            maxdelay = 16000;

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(1,N);
            
            Delay_time = zeros(1,N + maxdelay);
            Delay = 4000;


            for i =1:N

                if i == (N-1)/4
                    Delay = Delay * 2;
                end
                if i == (N-1)/2
                    Delay = Delay * 2;
                end
                if i == 3*(N-1)/4
                    Delay = Delay / 2;
                end

                Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
                if Delay ~= 0 && Delay_time(i) == 0 && i > 1
                    Delay_time(i) = Delay_time(i-1);
                end
%                 Delay_time(i) = Delay;
                muw = muw0*exp(-2*(Delay_time(i)/16000));
                MUW(i) = muw;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay,:) = Nabla1_store(i+Delay,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay,:) = Nabla2_store(i+Delay,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay,:) = Nabla3_store(i+Delay,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay,:) = Nabla4_store(i+Delay,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay,:) = Nabla5_store(i+Delay,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay,:) = Nabla6_store(i+Delay,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end

        %  delay; transmit gradient; 166 ; case 3; novss; sudden change delay
        function [e,Delay_time,MUW,obj] = DMANC_gradient_novss_case3(obj,xin,muw0)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            maxdelay = 16000;

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(1,N);
            
            Delay_time = zeros(1,N + maxdelay);
            Delay = 4000;


            for i =1:N

                if i == (N-1)/4
                    Delay = Delay * 2;
                end
                if i == (N-1)/2
                    Delay = Delay * 2;
                end
                if i == 3*(N-1)/4
                    Delay = Delay / 2;
                end

                Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
%                 Delay_time(i) = Delay;
                if Delay ~= 0 && Delay_time(i) == 0 && i > 1
                    Delay_time(i) = Delay_time(i-1);
                end
                muw = muw0;
                MUW(i) = muw;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay,:) = Nabla1_store(i+Delay,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay,:) = Nabla2_store(i+Delay,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay,:) = Nabla3_store(i+Delay,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay,:) = Nabla4_store(i+Delay,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay,:) = Nabla5_store(i+Delay,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay,:) = Nabla6_store(i+Delay,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end

        %  delay; transmit gradient; 166 ; case 5;
        %  each node has different delay
        function [e,Delay_time,Delay_store,MUW,obj] = DMANC_gradient_vss_case5(obj,xin,muw0)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            maxdelay = 16000;

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(6,N);
            
            Delay_time = zeros(6,N + maxdelay);
            Delay_store = zeros(6,N);
%             Delay1 = 4000;
%             Delay2 = 8000;
%             Delay3 = 4000;


            for i =1:N

%                 if i == (N-1)/4
%                     Delay1 = Delay1 * 2;  % 8000
%                     Delay2 = Delay2 / 2;  % 4000
%                     Delay3 = Delay3 * 4;  % 16000
%                 end
%                 if i == (N-1)/2
%                     Delay1 = Delay1 * 2;  % 16000
%                     Delay2 = Delay2 * 3;  % 12000
%                     Delay3 = Delay3 / 2;  % 8000
%                 end
%                 if i == 3*(N-1)/4
%                     Delay1 = Delay1 / 2;  % 8000
%                     Delay2 = Delay2 / 2;  % 6000
%                     Delay3 = Delay3 / 2;  % 4000
%                 end

%                 Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
%                 Delay_time(1,i+Delay1) = max(Delay_time(1,i+Delay1),Delay1);
%                 Delay_time(2,i+Delay2) = max(Delay_time(2,i+Delay2),Delay2);
%                 Delay_time(3,i+Delay3) = max(Delay_time(3,i+Delay3),Delay3);
%                 muw1 = muw0*exp(-2*(Delay_time(1,i)/16000));
%                 muw2 = muw0*exp(-2*(Delay_time(2,i)/16000));
%                 muw3 = muw0*exp(-2*(Delay_time(3,i)/16000));
%                 MUW(1,i) = muw1;
%                 MUW(2,i) = muw2;
%                 MUW(3,i) = muw3;

                Delay1 = round((sin(2*pi*0.05*1*(i/16000)-(pi/2))+1)*8000);
                Delay2 = round((sin(2*pi*0.05*2*(i/16000)-(pi/2))+1)*8000);
                Delay3 = round((sin(2*pi*0.05*3*(i/16000)-(pi/2))+1)*8000);
                Delay4 = round((sin(2*pi*0.05*4*(i/16000)-(pi/2))+1)*8000);
                Delay5 = round((sin(2*pi*0.05*5*(i/16000)-(pi/2))+1)*8000);
                Delay6 = round((sin(2*pi*0.05*6*(i/16000)-(pi/2))+1)*8000);

                Delay_store(1,i) = Delay1;
                Delay_store(2,i) = Delay2; 
                Delay_store(3,i) = Delay3; 
                Delay_store(4,i) = Delay4; 
                Delay_store(5,i) = Delay5; 
                Delay_store(6,i) = Delay6; 

                Delay_time(1,i+Delay1) = max(Delay_time(1,i+Delay1),Delay1);
                Delay_time(2,i+Delay2) = max(Delay_time(2,i+Delay2),Delay2);
                Delay_time(3,i+Delay3) = max(Delay_time(3,i+Delay3),Delay3);
                Delay_time(4,i+Delay4) = max(Delay_time(4,i+Delay4),Delay4);
                Delay_time(5,i+Delay5) = max(Delay_time(5,i+Delay5),Delay5);
                Delay_time(6,i+Delay6) = max(Delay_time(6,i+Delay6),Delay6);

                if Delay1 ~= 0 && Delay_time(1,i) == 0 && i > 1
                    Delay_time(1,i) = Delay_time(1,i-1);
                end
                if Delay2 ~= 0 && Delay_time(2,i) == 0 && i > 1
                    Delay_time(2,i) = Delay_time(2,i-1);
                end
                if Delay3 ~= 0 && Delay_time(3,i) == 0 && i > 1
                    Delay_time(3,i) = Delay_time(3,i-1);
                end
                if Delay4 ~= 0 && Delay_time(4,i) == 0 && i > 1
                    Delay_time(4,i) = Delay_time(4,i-1);
                end
                if Delay5 ~= 0 && Delay_time(5,i) == 0 && i > 1
                    Delay_time(5,i) = Delay_time(5,i-1);
                end
                if Delay6 ~= 0 && Delay_time(6,i) == 0 && i > 1
                    Delay_time(6,i) = Delay_time(6,i-1);
                end
                muw1 = muw0*exp(-2*(max([Delay_time(2,i) Delay_time(3,i) Delay_time(4,i) Delay_time(5,i) Delay_time(6,i)])/16000));
                muw2 = muw0*exp(-2*(max([Delay_time(1,i) Delay_time(3,i) Delay_time(4,i) Delay_time(5,i) Delay_time(6,i)])/16000));
                muw3 = muw0*exp(-2*(max([Delay_time(1,i) Delay_time(2,i) Delay_time(4,i) Delay_time(5,i) Delay_time(6,i)])/16000));
                muw4 = muw0*exp(-2*(max([Delay_time(1,i) Delay_time(2,i) Delay_time(3,i) Delay_time(5,i) Delay_time(6,i)])/16000));
                muw5 = muw0*exp(-2*(max([Delay_time(1,i) Delay_time(2,i) Delay_time(3,i) Delay_time(4,i) Delay_time(6,i)])/16000));
                muw6 = muw0*exp(-2*(max([Delay_time(1,i) Delay_time(2,i) Delay_time(3,i) Delay_time(4,i) Delay_time(5,i)])/16000));
                MUW(1,i) = muw1;
                MUW(2,i) = muw2;
                MUW(3,i) = muw3;
                MUW(4,i) = muw4;
                MUW(5,i) = muw5;
                MUW(6,i) = muw6;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay1,:) = Nabla1_store(i+Delay1,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay2,:) = Nabla2_store(i+Delay2,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay3,:) = Nabla3_store(i+Delay3,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay4,:) = Nabla4_store(i+Delay4,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay5,:) = Nabla5_store(i+Delay5,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay6,:) = Nabla6_store(i+Delay6,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw1*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw2*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw3*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw4*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw5*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw6*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end


        %  delay; transmit gradient; 166 ; case 5; sudden change delay,
        %  no vss
        %  each node has different delay
        function [e,Delay_time,MUW,obj] = DMANC_gradient_novss_case5(obj,xin,muw0,delaystore)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            maxdelay = 16000;

            Nabla1_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla2_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla3_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla4_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla5_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));
            Nabla6_store = zeros((N + maxdelay),(obj.wlen+obj.clen-1));

            Received_Nabla = zeros(obj.Numnode,(obj.wlen+obj.clen-1));

            MUW = zeros(3,N);
            
            Delay_time = zeros(6,N + maxdelay);
%             Delay1 = 4000;
%             Delay2 = 8000;
%             Delay3 = 4000;


            for i =1:N

%                 if i == (N-1)/4
%                     Delay1 = Delay1 * 2;  % 8000
%                     Delay2 = Delay2 / 2;  % 4000
%                     Delay3 = Delay3 * 4;  % 16000
%                 end
%                 if i == (N-1)/2
%                     Delay1 = Delay1 * 2;  % 16000
%                     Delay2 = Delay2 * 3;  % 12000
%                     Delay3 = Delay3 / 2;  % 8000
%                 end
%                 if i == 3*(N-1)/4
%                     Delay1 = Delay1 / 2;  % 8000
%                     Delay2 = Delay2 / 2;  % 6000
%                     Delay3 = Delay3 / 2;  % 4000
%                 end
% 
% %                 Delay_time(i+Delay) = max(Delay_time(i+Delay),Delay);
%                 Delay_time(1,i) = Delay1;
%                 Delay_time(2,i) = Delay2;
%                 Delay_time(3,i) = Delay3;
                
%                 Delay1 = round((sin(2*pi*0.1*(i/16000)-(pi/2))+1)*8000);
%                 Delay2 = round((sin(2*pi*0.2*(i/16000)-(pi/2))+1)*8000);
%                 Delay3 = round((sin(2*pi*0.3*(i/16000)-(pi/2))+1)*8000);
%                 Delay4 = round((sin(2*pi*0.4*(i/16000)-(pi/2))+1)*8000);
%                 Delay5 = round((sin(2*pi*0.5*(i/16000)-(pi/2))+1)*8000);
%                 Delay6 = round((sin(2*pi*0.6*(i/16000)-(pi/2))+1)*8000);

                Delay1 = delaystore(1,i);
                Delay2 = delaystore(2,i);
                Delay3 = delaystore(3,i);
                Delay4 = delaystore(4,i);
                Delay5 = delaystore(5,i);
                Delay6 = delaystore(6,i);

                Delay_time(1,i+Delay1) = max(Delay_time(1,i+Delay1),Delay1);
                Delay_time(2,i+Delay2) = max(Delay_time(2,i+Delay2),Delay2);
                Delay_time(3,i+Delay3) = max(Delay_time(3,i+Delay3),Delay3);
                Delay_time(4,i+Delay4) = max(Delay_time(4,i+Delay4),Delay4);
                Delay_time(5,i+Delay5) = max(Delay_time(5,i+Delay5),Delay5);
                Delay_time(6,i+Delay6) = max(Delay_time(6,i+Delay6),Delay6);

                if Delay1 ~= 0 && Delay_time(1,i) == 0 && i > 1
                    Delay_time(1,i) = Delay_time(1,i-1);
                end
                if Delay2 ~= 0 && Delay_time(2,i) == 0 && i > 1
                    Delay_time(2,i) = Delay_time(2,i-1);
                end
                if Delay3 ~= 0 && Delay_time(3,i) == 0 && i > 1
                    Delay_time(3,i) = Delay_time(3,i-1);
                end
                if Delay4 ~= 0 && Delay_time(4,i) == 0 && i > 1
                    Delay_time(4,i) = Delay_time(4,i-1);
                end
                if Delay5 ~= 0 && Delay_time(5,i) == 0 && i > 1
                    Delay_time(5,i) = Delay_time(5,i-1);
                end
                if Delay6 ~= 0 && Delay_time(6,i) == 0 && i > 1
                    Delay_time(6,i) = Delay_time(6,i-1);
                end
                
                muw1 = muw0;
                muw2 = muw0;
                muw3 = muw0;
                MUW(1,i) = muw1;
                MUW(2,i) = muw2;
                MUW(3,i) = muw3;
                
                xc = [xin(i) xc(1:(end-1))];   % update reference vector

                Nabla1_store(i+Delay1,:) = Nabla1_store(i+Delay1,:) + obj.Nabla(1,:);
                Nabla2_store(i+Delay2,:) = Nabla2_store(i+Delay2,:) + obj.Nabla(2,:);
                Nabla3_store(i+Delay3,:) = Nabla3_store(i+Delay3,:) + obj.Nabla(3,:);
                Nabla4_store(i+Delay4,:) = Nabla4_store(i+Delay4,:) + obj.Nabla(4,:);
                Nabla5_store(i+Delay5,:) = Nabla5_store(i+Delay5,:) + obj.Nabla(5,:);
                Nabla6_store(i+Delay6,:) = Nabla6_store(i+Delay6,:) + obj.Nabla(6,:);


                Received_Nabla(1,:) = Nabla1_store(i,:);
                Received_Nabla(2,:) = Nabla2_store(i,:);
                Received_Nabla(3,:) = Nabla3_store(i,:);
                Received_Nabla(4,:) = Nabla4_store(i,:);
                Received_Nabla(5,:) = Nabla5_store(i,:);
                Received_Nabla(6,:) = Nabla6_store(i,:);

                
                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + muw0*(obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + muw0*(obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + muw0*(obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + muw0*(obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,Received_Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + muw0*(obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,Received_Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,Received_Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,Received_Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,Received_Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,Received_Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + muw0*(obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end)+ ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end));


                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = xf(6,:)*e(6,i);


            end
        end

       

    end
end