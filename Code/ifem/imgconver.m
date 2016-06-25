function imconver
  [a,result] = system('find . -name \*.fig -print');
  result = strsplit(result, '\n');
  for i=1:length(result)
    file = result{i}
    if strcmp(file, '')
      continue
    end
    figure('visible' ,'off')
    fig = openfig(file, 'invisible');
    export_fig(fig, strrep(file, '.fig' ,'_ef'),'-pdf', '-png', '-eps', '-painters');
    %saveas(fig, strrep(file, '.fig', '.eps'))
  end
end
