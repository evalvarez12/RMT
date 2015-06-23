#RMT in julia
#Eduardo Villaseñor - evalvarez12@gmail.com

#   itpp::Mat<std::complex<double> > RandomGUE(int const dim, std::string normalization="sigma_offdiag=1", double const percentage_away=0.1){ //{{{
#     if (normalization=="sigma_offdiag=1"){
#       return RandomGUEDeltaOne(dim);
#     }	
#     else if (normalization=="unfolded mean_level_spacing=1"){
#       itpp::Mat<std::complex<double> > U(dim, dim), tmp(dim,dim);
#       itpp::Vec<std::complex<double> > vec1(dim);
#       itpp::Vec<double> eigenvalues(dim);
#       FlatSpectrumGUE(U, eigenvalues);
#       for (int i=0; i<dim; i++){
#         vec1=itpp::elem_mult(conj(U.get_col(i)), to_cvec(eigenvalues));
#         for (int j=i; j<dim; j++){
#           tmp(i,j)=vec1*U.get_col(j);
#           if (i<j){tmp(j,i)=conj(tmp(i,j));}
#         }
#       }
#       return tmp;
#       
#       
#       
#         itpp::Mat<std::complex<double> > RandomGUEDeltaOne(int const dim){ //{{{
#     itpp::Mat<std::complex<double> > temp(dim, dim);
#     temp=itpp::randn_c(dim,dim);
#     return sqrt(1/2.)*(temp+temp.hermitian_transpose());
#   }

function RandomGUEDeltaOne(dim) 
  A= random(dim,dim)
  