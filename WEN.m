function [Envelope] = WEN(m,width)
FW_Envelope=slider(m,width);
diff=length(m)-length(FW_Envelope);
back_m=m(end:-1:1);
BW_Envelope=slider(back_m,width);
BW_Envelope=BW_Envelope(end:-1:1);
% FW_Envelope=[FW_Envelope,BW_Envelope(end-diff+1:end)];
% BW_Envelope=[FW_Envelope(1:diff),BW_Envelope];
FW_Envelope=[FW_Envelope,nan(1,diff)];
BW_Envelope=[nan(1,diff),BW_Envelope];
Envelope = zeros(size(BW_Envelope));
for i=1:length(FW_Envelope)
    Envelope(i)=min(FW_Envelope(i),BW_Envelope(i));
end
end