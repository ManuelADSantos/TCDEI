function [result_rms,wrong_type]=RecursiveTest(EstFunc,NumberIter,NoiseStdValue,ConicType,M);

Cut=0.1;
N_func=max(size(EstFunc));
NumberIter;
N=max(size(NoiseStdValue));
result_rms=zeros(N_func,N);
wrong_type=result_rms;
for i=1:1:N
 error_vec=zeros(N_func,NumberIter);
 acc=zeros(N_func,1);
 for j=1:1:NumberIter
  [Points,gt,dummy]=CurveGenerator(NoiseStdValue(i),M,ConicType);
  for k=1:1:N_func
   omega=feval(EstFunc{k},Points);
   if sign(omega(1,2)^2-omega(1,1)*omega(2,2)) ~= sign(gt(1,2)^2-gt(1,1)*gt(2,2))
    acc(k)=acc(k)+1;
   end;
   error_vec(k,j)=ConicCompare(gt,omega);
  end;
 end;
 %Eliminate too bad results
 error_vec=sort(error_vec,2);
 error_vec=error_vec(:,round(Cut*NumberIter+1):NumberIter-round(Cut*NumberIter));
 result_rms(:,i)=sqrt(mean(error_vec.^2,2)/NumberIter);
 wrong_type(:,i)=acc./NumberIter;
end;

  
