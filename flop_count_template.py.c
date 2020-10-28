#ifdef PYTHON
#define ONELINE_FUNCDEF(name, body) def name(): return (body)
#else
#define ONELINE_FUNCDEF(name, body) static float name(){ return (body);}
#endif


// see Dphi_flt.c

#ifdef REPR_ADJOINT
ONELINE_FUNCDEF(matrix_mul_flops, 4*NF*NF - 2*NF);
#else
ONELINE_FUNCDEF(matrix_mul_flops, 8*NF*NF - 2*NF);
#endif

ONELINE_FUNCDEF(real_op_flop , 1);
ONELINE_FUNCDEF(complex_sum_flops , 2);
ONELINE_FUNCDEF(complex_mul_flops , 6);
ONELINE_FUNCDEF(complex_mulR_flops , 2); // capital R for legibility (hopefully);
ONELINE_FUNCDEF(complex_norm_flops , 3); // only re ;

// single vector
ONELINE_FUNCDEF(vector_sum_flops , complex_sum_flops()*NF);
ONELINE_FUNCDEF(vector_mul_flops , complex_mulR_flops()*NF); // mul uses mulR ;
ONELINE_FUNCDEF(vector_norm_flops , complex_norm_flops()*NF+real_op_flop()*NF); // + NF sums (not NF-1)

// single spinor
ONELINE_FUNCDEF(spinor_norm_flops , 4*vector_norm_flops() );
ONELINE_FUNCDEF(spinor_sum_flops , 4*vector_sum_flops());
ONELINE_FUNCDEF(spinor_mul_flops , 4*vector_mul_flops()); // uses mulR
ONELINE_FUNCDEF(spinor_mul_add_assign_flops , spinor_mul_flops() + spinor_sum_flops());

// spinor field related flop count
ONELINE_FUNCDEF(site_spinor_field_sqnorm_f_flops , 1.0/2 * spinor_norm_flops()); // mult + add
ONELINE_FUNCDEF(site_spinor_field_mul_add_assign_f_flops , 1.0/2 * spinor_mul_add_assign_flops());
ONELINE_FUNCDEF(site_spinor_field_mul_f_flops , 1.0/2 * spinor_mul_flops());
ONELINE_FUNCDEF(site_spinor_field_sub_f_flops , 1.0/2 * spinor_sum_flops());
ONELINE_FUNCDEF(site_spinor_field_prod_re_flop , site_spinor_field_sqnorm_f_flops());
ONELINE_FUNCDEF(site_spinor_field_add_assign_f_flops , 1.0/2 * spinor_sum_flops());
ONELINE_FUNCDEF(site_spinor_field_minus_f_flops , 1.0/2 * spinor_mul_flops()); // as a mul

// see Dphi.c
ONELINE_FUNCDEF(site_Dphi_flops , 1.0/2*(16*matrix_mul_flops() +   \
                                          45*vector_sum_flops()));

ONELINE_FUNCDEF(site_Cphi_flops , 1.0/2*( 8*matrix_mul_flops() +   \
                                           4*vector_sum_flops() +   \
                                           spinor_mul_add_assign_flops()));

ONELINE_FUNCDEF(site_Cphi_assign_flops , 1.0/2*( 8*matrix_mul_flops() + \
                                                  4*vector_sum_flops() + \
                                                  spinor_mul_add_assign_flops() + \
                                                  spinor_sum_flops()));



ONELINE_FUNCDEF( N , 2*NF);
ONELINE_FUNCDEF( forward_substitution_loop_flop , N()*(N()-1)/2 * ((3 * real_op_flop()) + // n = i*(i+1)/2+i) \
                                                     // _complex_mul_sub_assign(); )
                                                     (complex_mul_flops() + complex_sum_flops()) + \
                                                     // _complex_mul_sub_assign();)
                                                     (complex_mul_flops() + complex_sum_flops())));

#ifdef PYTHON
def backward_substitution_loop_flop():
    res = 0
#else
float backward_substitution_loop_flop(){
    int res = 0;
#endif

#ifdef PYTHON
    for i in range(N()-1,-1,-1):
#else
    for(int i = N()-1; i >= 0; i--){
#endif
        res += 3; // n = i*(i+1)/2+i);
        res += complex_mulR_flops();//_complex_mulr());
        res += complex_mulR_flops();//_complex_mulr());
#ifdef PYTHON
        for k in range(i+1,N()):
#else
        for(int k = i+1; k < N(); k++){
#endif
            res += 3; // n = k*(k+1)/2+i);

            // _complex_mul_sub_assign());
            res += complex_mul_flops() + complex_sum_flops();
            // _complex_mul_sub_assign());
            res += complex_mul_flops() + complex_sum_flops();

#ifndef PYTHON
        }
    }
#endif
    return res;
#ifndef PYTHON
}
#endif

ONELINE_FUNCDEF( site_Cphi_inv_flops , 1.0/2*(forward_substitution_loop_flop() + \
                                      backward_substitution_loop_flop()));

ONELINE_FUNCDEF( site_Cphi_inv_assign_flops , 1.0/2*(forward_substitution_loop_flop() + \
                                             backward_substitution_loop_flop() + \
                                             spinor_sum_flops()));

#ifdef PYTHON
def site_g5Cphi_eopre_sq_flops():
#else
float site_g5Cphi_eopre_sq_flops(){
    float site_g5Cphi_eopre_flops;
#endif
    site_g5Cphi_eopre_flops = (site_Dphi_flops() +
                               site_Cphi_inv_flops() +
                               site_Dphi_flops() +
                               site_spinor_field_minus_f_flops() +
                               site_Cphi_assign_flops());

    return 2*site_g5Cphi_eopre_flops;

#ifndef PYTHON
}
#endif

#ifdef PYTHON
def cg_out_of_loop_flops_per_site(site_operator_flops):
#else
float cg_out_of_loop_flops_per_site(float site_operator_flops){
#endif
    return (site_spinor_field_sqnorm_f_flops() +
           site_operator_flops +
           site_spinor_field_mul_add_assign_f_flops() +
           site_spinor_field_sub_f_flops() +
           site_spinor_field_sqnorm_f_flops() +
           site_operator_flops +
           site_spinor_field_mul_add_assign_f_flops+
           site_spinor_field_sub_f_flops() +
           2*site_spinor_field_sqnorm_f_flops());
#ifndef PYTHON
}
#endif

#ifdef PYTHON
def cg_iteration_flops_per_site(site_operator_flops):
#else
float cg_iteration_flops_per_site(float site_operator_flops) {
#endif
    return (site_operator_flops +
            site_spinor_field_prod_re_flop() +           
            site_spinor_field_mul_add_assign_f_flops() + 
            site_spinor_field_mul_add_assign_f_flops() + 
            site_spinor_field_sqnorm_f_flops() +         
            site_spinor_field_mul_f_flops() +            
            site_spinor_field_add_assign_f_flops() +     
            site_spinor_field_mul_f_flops() +            
            site_spinor_field_add_assign_f_flops());
#ifndef PYTHON
}
#endif
