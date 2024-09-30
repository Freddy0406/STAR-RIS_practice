function [loc,reflect_element_num,element_width,element_size]=configRIS(wave_len)

    loc = [50,10];
    reflect_element_num = 64;       % Must be perfect square
    element_width = wave_len/4;
    element_size = element_width^2;

end