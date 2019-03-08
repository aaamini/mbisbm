function [ yh ] = soft_to_hard( yt )
%SOFTTOHARD Summary of this function goes here
%   Detailed explanation goes here

yh = label_vec2mat(label_mat2vec(yt));



