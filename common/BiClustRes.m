classdef BiClustRes
   properties
      tau, lab, cm, nmi, miss_rate
      name 
   end
   methods
       function ob = BiClustRes(data, name)
           if nargin < 2
               ob.name = 'generic';
           else
               ob.name = name;
           end
           d1 = size(data{1},2);
           d2 = size(data{2},2);
           if d1 ~= d2, error('labels should be of the same type.'), end
           
           if d1 == 1
               ob.lab = {};
               ob.lab{1} = data{1};
               ob.lab{2} = data{2};
               ob.tau = {};
               ob.tau{1} = label_vec2mat(ob.lab{1});
               ob.tau{2} = label_vec2mat(ob.lab{2});
           else
               ob.tau = {};
               ob.tau{1} = data{1};
               ob.tau{2} = data{2};
               ob.lab = {};
               ob.lab{1} = label_mat2vec(ob.tau{1});
               ob.lab{2} = label_mat2vec(ob.tau{2});
           end
                      
           
           tru_lab = BiClustRes.setget_tru_lab;
           ob.cm = confusionmat(            [tru_lab{1}(:); tru_lab{2}(:)], [ob.lab{1}; ob.lab{2}] );
           ob.nmi = compute_mutual_info(    [tru_lab{1}(:); tru_lab{2}(:)], [ob.lab{1}; ob.lab{2}] );
           
           if BiClustRes.setget_comp_miss
                ob.miss_rate = misclass_rate(ob.cm)*100;
           end
       end
       
       function disp(ob)
           fprintf('%-30s %5.2f    %5.2f%%\n', ob.name, ob.nmi, ob.miss_rate)
       end
       
       function out = labels_cat(ob)
           out = [ob.lab{1}; ob.lab{2}];
       end
   end
   methods (Static)
      function out = setget_tru_lab(data)
         persistent tru_lab;
         if nargin
            tru_lab = data;
         end
         out = tru_lab;
      end
      
       function out = setget_comp_miss(data)
         persistent comp_miss;
         if nargin
            comp_miss = data;
         end
         out = comp_miss;
      end
   end
end