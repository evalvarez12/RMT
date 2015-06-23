#RMT in julia
#Eduardo Villaseñor - evalvarez12@gmail.com


function RandomGUE(n::Integer) 
  A = randn(n, n) + im*randn(n, n)
  normalization = sqrt(4*n)
  return (A + A') / normalization
end
  
 