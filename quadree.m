function quadtree(file_path)
    % Wczytaj plik
    f = imread(file_path);

    % Wybieramy liczbę p = nextpow2(max(size(f)) tak aby 2^p > max(size(f)
    q = 2^nextpow2(max(size(f)));

    % Rozmiary obrazu
    [m, n]=size(f);

    % Dodanie ramki w okół obrazu
    f = padarray(f,[q-m,q-n],'post');

    % Minimalny podział (musi być potęgą dwójki)
    mindim = 2;

    % Dzieli kwadrat na cztery mniejsze kwadraty i testuje czy każdy z nich spełnia kryteria
    % Jeśli kwadrat spełnia kryteria nie jest już dzielony, jeśli nie podział następuje dalej
    s = qtdecomp(f, @split, mindim, @predicate);

    % Konwersja sparse macierzy na full macierz
    lmax = full(max(s(:)));

    % Inicjalizacji obrazu wyjściowego zerami
    g = zeros(size(f));

    % Inicjalizacja markerów zerami
    marker = zeros(size(f));

    for k=1:lmax
        [vals, r, c] = qtgetblk(f,s,k);
        if ~isempty(vals)
            for i=1:length(r)
                xlow = r(i);ylow=c(i);
                xhigh = xlow+k-1;
                yhigh = ylow+k-1;
                region = f(xlow:xhigh,ylow:yhigh);
                flag = feval(@predicate,region);
                if flag
                     g(xlow:xhigh,ylow:yhigh) = 1;
                     marker(xlow,ylow) = 1;
                end
            end
        end
    end
    g = bwlabel(imreconstruct(marker, g));

    % Przywrócenie wymiarów wejścowych
    g = g(1:m,1:n);
    f = f(1:m,1:n);

    figure, imshow(f),title('Original Image');
    figure, imshow(g),title('Segmented Image');
end
