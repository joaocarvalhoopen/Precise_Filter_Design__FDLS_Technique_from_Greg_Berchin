#ifndef CUSTOM_FILTER_H_
#define CUSTOM_FILTER_H_

//generated_C_.h_source_code_filter_file


float const a_1 = -1.844215;
float const a_2 = 0.985514;

float const b_0 = 0.303389;
float const b_1 = -0.576379;
float const b_2 = 0.303625;

float p_0 = (float) 0.0;
float p_1 = (float) 0.0;
float p_2 = (float) 0.0;

void custom_filter_init(){
    float p_0 = 0.0;
    float p_1 = 0.0;
    float p_2 = 0.0;
}

float custom_c_filter(float inputValue){
    p_0 = inputValue  + (-a_1*p_1) + (-a_2*p_2);
    float outputValue =  + (b_0*p_0) + (b_1*p_1) + (b_2*p_2);
    p_2 = p_1;
    p_1 = p_0;
    return outputValue;
}

#endif
