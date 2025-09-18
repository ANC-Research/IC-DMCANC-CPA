%% Conventional Multichannel ANC (1 x K x M)  SIMO
% 1 reference, K secondary sources, M error sensors


classdef McANC_FxLMS_SIMO
    properties
        Wc    % control filter (K * wlen)
        wlen  % length of control filter
        SecP  % secondary path estimates (M * k * slen)
        slen  % length of secondary path
        Xd    % Disturbance
        yc    % control filter output
        Nums  % number of secondary sources, K
        Nume  % number of errors, M

    end

    methods
        % initialize (N:duration; PrimaryPath: M x plen)
        function obj = McANC_FxLMS_SIMO(wLen,SecondaryPath,sLen,Nums,Nume,N,Dis)
            obj.wlen = wLen;
            obj.Wc   = zeros(Nums,wLen);
            obj.SecP = SecondaryPath;
            obj.slen = sLen;
            obj.yc   = zeros(Nums,N);
%             obj.Xd   = zeros(Nume,N);
            obj.Xd   = Dis;
            obj.Nume = Nume;
            obj.Nums = Nums;
%             for i = 1:Nume
%                 obj.Xd(i,:) = filter(PrimaryPath(i,:),1,xin);  % disturbance received by error sensors
%             end
        end
        % 1 x 2 x 2
        function [e,obj] = McFxLMS_SIMO_122(obj,xin,muw)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:));      % y1' received by error 1
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:));      % y2' received by error 2
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                
                xs = [xin(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i));  % update control filter 1
                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i));  % update control filter 2
            end
        end

                % 1 x 4 x 4
        function [e,obj] = McFxLMS_SIMO_144(obj,xin,muw)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4
                
                xs = [xin(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41

                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42

                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43

                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i));  % update control filter 1
                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i));  % update control filter 2
                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i));  % update control filter 1
                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i));  % update control filter 2
            end
        end
            

        function [e,obj] = McFxLMS_SIMO_144_initialCF(obj,xin,muw,CF)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            obj.Wc = CF;

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update


                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4

                xs = [xin(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41


                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42


                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43


                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44


                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i));  % update control filter 1

                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i));  % update control filter 2

                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i));  % update control filter 3

                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i));  % update control filter 4


            end
        end
        
        function [e,obj] = McFxLMS_SIMO_144_fixed(obj,xin,CF)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter

            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path
            obj.Wc = CF;

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4


                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update


                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:));

                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4
                

            end
        end


        function [e,obj] = McFxLMS_SIMO_166(obj,xin,muw)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);         % control output y5
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);         % control output y6

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));

                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));

                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4
                e(5,i) = obj.Xd(5,i)-y5;                    % error 5
                e(6,i) = obj.Xd(6,i)-y6;                    % error 6
                
                xs = [xin(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41
                xf(1,5,:) = [sum(xs.*reshape(obj.SecP(5,1,:),[1,obj.slen])) reshape(xf(1,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s51
                xf(1,6,:) = [sum(xs.*reshape(obj.SecP(6,1,:),[1,obj.slen])) reshape(xf(1,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s61

                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42
                xf(2,5,:) = [sum(xs.*reshape(obj.SecP(5,2,:),[1,obj.slen])) reshape(xf(2,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s52
                xf(2,6,:) = [sum(xs.*reshape(obj.SecP(6,2,:),[1,obj.slen])) reshape(xf(2,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s62

                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43
                xf(3,5,:) = [sum(xs.*reshape(obj.SecP(5,3,:),[1,obj.slen])) reshape(xf(3,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s53
                xf(3,6,:) = [sum(xs.*reshape(obj.SecP(6,3,:),[1,obj.slen])) reshape(xf(3,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s63

                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44
                xf(4,5,:) = [sum(xs.*reshape(obj.SecP(5,4,:),[1,obj.slen])) reshape(xf(4,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s54
                xf(4,6,:) = [sum(xs.*reshape(obj.SecP(6,4,:),[1,obj.slen])) reshape(xf(4,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s64

                xf(5,1,:) = [sum(xs.*reshape(obj.SecP(1,5,:),[1,obj.slen])) reshape(xf(5,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s15
                xf(5,2,:) = [sum(xs.*reshape(obj.SecP(2,5,:),[1,obj.slen])) reshape(xf(5,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s25
                xf(5,3,:) = [sum(xs.*reshape(obj.SecP(3,5,:),[1,obj.slen])) reshape(xf(5,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s35
                xf(5,4,:) = [sum(xs.*reshape(obj.SecP(4,5,:),[1,obj.slen])) reshape(xf(5,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s45
                xf(5,5,:) = [sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen])) reshape(xf(5,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s55
                xf(5,6,:) = [sum(xs.*reshape(obj.SecP(6,5,:),[1,obj.slen])) reshape(xf(5,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s65

                xf(6,1,:) = [sum(xs.*reshape(obj.SecP(1,6,:),[1,obj.slen])) reshape(xf(6,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s16
                xf(6,2,:) = [sum(xs.*reshape(obj.SecP(2,6,:),[1,obj.slen])) reshape(xf(6,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s26
                xf(6,3,:) = [sum(xs.*reshape(obj.SecP(3,6,:),[1,obj.slen])) reshape(xf(6,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s36
                xf(6,4,:) = [sum(xs.*reshape(obj.SecP(4,6,:),[1,obj.slen])) reshape(xf(6,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s46
                xf(6,5,:) = [sum(xs.*reshape(obj.SecP(5,6,:),[1,obj.slen])) reshape(xf(6,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s56
                xf(6,6,:) = [sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen])) reshape(xf(6,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s66

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(1,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(1,6,:),[1,obj.wlen])*e(6,i));  % update control filter 1

                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(2,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(2,6,:),[1,obj.wlen])*e(6,i));  % update control filter 2

                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(3,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(3,6,:),[1,obj.wlen])*e(6,i));  % update control filter 3

                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(4,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(4,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4

                obj.Wc(5,:) = obj.Wc(5,:) + muw*(reshape(xf(5,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(5,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(5,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(5,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(5,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(5,6,:),[1,obj.wlen])*e(6,i));  % update control filter 5

                obj.Wc(6,:) = obj.Wc(6,:) + muw*(reshape(xf(6,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(6,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(6,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(6,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(6,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(6,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4
            end
        end
    
        function [e,obj] = McFxLMS_SIMO_166_initialCF(obj,xin,muw,CF)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter
            xs = zeros(1,obj.slen);                         % x buffer for filter secondary path
            xf = zeros(obj.Nums,obj.Nume,obj.wlen);         % filtered x buffer for update
            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path

            obj.Wc = CF;

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);         % control output y5
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);         % control output y6

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));

                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));

                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4
                e(5,i) = obj.Xd(5,i)-y5;                    % error 5
                e(6,i) = obj.Xd(6,i)-y6;                    % error 6
                
                xs = [xin(i) xs(1:(end-1))];                                        % update filter x buffer
                xf(1,1,:) = [sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen])) reshape(xf(1,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s11
                xf(1,2,:) = [sum(xs.*reshape(obj.SecP(2,1,:),[1,obj.slen])) reshape(xf(1,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s21
                xf(1,3,:) = [sum(xs.*reshape(obj.SecP(3,1,:),[1,obj.slen])) reshape(xf(1,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s31
                xf(1,4,:) = [sum(xs.*reshape(obj.SecP(4,1,:),[1,obj.slen])) reshape(xf(1,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s41
                xf(1,5,:) = [sum(xs.*reshape(obj.SecP(5,1,:),[1,obj.slen])) reshape(xf(1,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s51
                xf(1,6,:) = [sum(xs.*reshape(obj.SecP(6,1,:),[1,obj.slen])) reshape(xf(1,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s61

                xf(2,1,:) = [sum(xs.*reshape(obj.SecP(1,2,:),[1,obj.slen])) reshape(xf(2,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s12
                xf(2,2,:) = [sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen])) reshape(xf(2,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s22
                xf(2,3,:) = [sum(xs.*reshape(obj.SecP(3,2,:),[1,obj.slen])) reshape(xf(2,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s32
                xf(2,4,:) = [sum(xs.*reshape(obj.SecP(4,2,:),[1,obj.slen])) reshape(xf(2,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s42
                xf(2,5,:) = [sum(xs.*reshape(obj.SecP(5,2,:),[1,obj.slen])) reshape(xf(2,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s52
                xf(2,6,:) = [sum(xs.*reshape(obj.SecP(6,2,:),[1,obj.slen])) reshape(xf(2,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s62

                xf(3,1,:) = [sum(xs.*reshape(obj.SecP(1,3,:),[1,obj.slen])) reshape(xf(3,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s13
                xf(3,2,:) = [sum(xs.*reshape(obj.SecP(2,3,:),[1,obj.slen])) reshape(xf(3,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s23
                xf(3,3,:) = [sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen])) reshape(xf(3,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s33
                xf(3,4,:) = [sum(xs.*reshape(obj.SecP(4,3,:),[1,obj.slen])) reshape(xf(3,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s43
                xf(3,5,:) = [sum(xs.*reshape(obj.SecP(5,3,:),[1,obj.slen])) reshape(xf(3,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s53
                xf(3,6,:) = [sum(xs.*reshape(obj.SecP(6,3,:),[1,obj.slen])) reshape(xf(3,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s63

                xf(4,1,:) = [sum(xs.*reshape(obj.SecP(1,4,:),[1,obj.slen])) reshape(xf(4,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s14
                xf(4,2,:) = [sum(xs.*reshape(obj.SecP(2,4,:),[1,obj.slen])) reshape(xf(4,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s24
                xf(4,3,:) = [sum(xs.*reshape(obj.SecP(3,4,:),[1,obj.slen])) reshape(xf(4,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s34
                xf(4,4,:) = [sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen])) reshape(xf(4,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s44
                xf(4,5,:) = [sum(xs.*reshape(obj.SecP(5,4,:),[1,obj.slen])) reshape(xf(4,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s54
                xf(4,6,:) = [sum(xs.*reshape(obj.SecP(6,4,:),[1,obj.slen])) reshape(xf(4,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s64

                xf(5,1,:) = [sum(xs.*reshape(obj.SecP(1,5,:),[1,obj.slen])) reshape(xf(5,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s15
                xf(5,2,:) = [sum(xs.*reshape(obj.SecP(2,5,:),[1,obj.slen])) reshape(xf(5,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s25
                xf(5,3,:) = [sum(xs.*reshape(obj.SecP(3,5,:),[1,obj.slen])) reshape(xf(5,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s35
                xf(5,4,:) = [sum(xs.*reshape(obj.SecP(4,5,:),[1,obj.slen])) reshape(xf(5,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s45
                xf(5,5,:) = [sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen])) reshape(xf(5,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s55
                xf(5,6,:) = [sum(xs.*reshape(obj.SecP(6,5,:),[1,obj.slen])) reshape(xf(5,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s65

                xf(6,1,:) = [sum(xs.*reshape(obj.SecP(1,6,:),[1,obj.slen])) reshape(xf(6,1,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s16
                xf(6,2,:) = [sum(xs.*reshape(obj.SecP(2,6,:),[1,obj.slen])) reshape(xf(6,2,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s26
                xf(6,3,:) = [sum(xs.*reshape(obj.SecP(3,6,:),[1,obj.slen])) reshape(xf(6,3,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s36
                xf(6,4,:) = [sum(xs.*reshape(obj.SecP(4,6,:),[1,obj.slen])) reshape(xf(6,4,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s46
                xf(6,5,:) = [sum(xs.*reshape(obj.SecP(5,6,:),[1,obj.slen])) reshape(xf(6,5,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s56
                xf(6,6,:) = [sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen])) reshape(xf(6,6,1:(obj.wlen-1)),[1,obj.wlen-1])];      % filtered x by s66

                obj.Wc(1,:) = obj.Wc(1,:) + muw*(reshape(xf(1,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(1,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(1,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(1,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(1,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(1,6,:),[1,obj.wlen])*e(6,i));  % update control filter 1

                obj.Wc(2,:) = obj.Wc(2,:) + muw*(reshape(xf(2,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(2,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(2,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(2,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(2,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(2,6,:),[1,obj.wlen])*e(6,i));  % update control filter 2

                obj.Wc(3,:) = obj.Wc(3,:) + muw*(reshape(xf(3,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(3,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(3,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(3,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(3,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(3,6,:),[1,obj.wlen])*e(6,i));  % update control filter 3

                obj.Wc(4,:) = obj.Wc(4,:) + muw*(reshape(xf(4,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(4,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(4,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(4,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(4,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(4,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4

                obj.Wc(5,:) = obj.Wc(5,:) + muw*(reshape(xf(5,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(5,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(5,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(5,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(5,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(5,6,:),[1,obj.wlen])*e(6,i));  % update control filter 5

                obj.Wc(6,:) = obj.Wc(6,:) + muw*(reshape(xf(6,1,:),[1,obj.wlen])*e(1,i) + reshape(xf(6,2,:),[1,obj.wlen])*e(2,i) + ...
                                                 reshape(xf(6,3,:),[1,obj.wlen])*e(3,i) + reshape(xf(6,4,:),[1,obj.wlen])*e(4,i) + ...
                                                 reshape(xf(6,5,:),[1,obj.wlen])*e(5,i) + reshape(xf(6,6,:),[1,obj.wlen])*e(6,i));  % update control filter 4
            end
        end
        
        function [e,obj] = McFxLMS_SIMO_166_fixed(obj,xin,CF)
            N  = length(xin);                               % duration
            e  = zeros(obj.Nume,N);                         % error signal M x duration
            xc = zeros(1,obj.wlen);                         % x buffer for control filter

            ys = zeros(obj.Nums,obj.slen);                  % y buffer for filter secondary path
            obj.Wc = CF;

            for i = 1:N
                xc = [xin(i) xc(1:(end-1))];                % update control filter x buffer
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);         % control output y1
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);         % control output y2
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);         % control output y3
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);         % control output y4
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);         % control output y5
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);         % control output y6

                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % y received by error
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));

                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));

                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));

                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));

                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));

                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                
                e(1,i) = obj.Xd(1,i)-y1;                    % error 1
                e(2,i) = obj.Xd(2,i)-y2;                    % error 2
                e(3,i) = obj.Xd(3,i)-y3;                    % error 3
                e(4,i) = obj.Xd(4,i)-y4;                    % error 4
                e(5,i) = obj.Xd(5,i)-y5;                    % error 5
                e(6,i) = obj.Xd(6,i)-y6;                    % error 6
                

            end
        end

    end
end


