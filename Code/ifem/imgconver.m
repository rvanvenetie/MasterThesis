function imconver
  [a,result] = system('find . -name \*.fig -print');
  result = strsplit(result, '\n');
  for i=1:length(result)
    file = result{i}
    if strcmp(file, '')
      continue
    end
    % Delete all bullshit files
    delete(strrep(file, '.fig', '.png'), strrep(file, '.fig', '.eps'), strrep(file, '.fig', '.pdf'));
    delete(strrep(file, '.fig', '_ef.png'), strrep(file, '.fig', '_ef.eps'), strrep(file, '.fig', '_ef.pdf'));
    warning('off','last');

    % Do not open this figure
    figure('visible' ,'off')
    fig = openfig(file, 'invisible');
    set(fig, 'Color', 'none'); % Set background transparent
    export_fig(fig, strrep(file, '.fig' ,''),'-pdf', '-painters'); % Export to PDF using export_fig
    %saveas(fig, strrep(file, '.fig', '.eps'))
  end
end
