%% Limpar os Códigos Antigos
% Aqui é realizado a limpeza do código assim como das imagens ante-
% riores.
clc
clear all
close all

%% Carregar a nova Imagem

load('Mama1.mat'); % Aqui é realizado o carregamento da imagem.
Mama=Mama1;

[l c] = size(Mama); %Determinação do tamanho da imagem original.
% imagesc(Mama)

BW = edge(Mama,'Sobel'); %Detecção das bordas pelo metodo de Sobel.

%% Criar histograma
% A criação de histogramas e importante é importante para se saber que
% existem duas temporaturas diferentes para o plano de fundo assim como
% saber qual o melhor tresholding para se distinguir a imagem de fundo do
% corpo da pessoa.

%Criação de um histograma com divisões de 0.01, isto é o limite da camera.
figure
A = histogram(Mama,round((max(max(Mama))-min(min(Mama)))*1000));
Hist = [A.BinCounts; A.BinEdges(2:end)];
Iblur = imgaussfilt(Hist(1,:), 100);
HistNew = [Iblur; A.BinEdges(2:end)];

%Tentativa para a detecção dos picos de modo a se separar imagem da pessoa
%e a imagem de fundo.

% [B C] = findpeaks(Hist(1,:),'NPeaks', 10);
% [B C] = findpeaks(Iblur(1,:),Hist(2,:),'NPeaks', 5,'MinPeakWidth',0.10 );
[XX YY] = findpeaks(Iblur(1,:),Hist(2,:),'NPeaks', 5,'MinPeakProminence',1 );
figure
plot(HistNew(2,:),HistNew(1,:));

% J = imwarp(I,tform);
% tform = affine2d([1 0 0; .5 1 0; 0 0 1]);

%% Encontrar os pontos de simetria da imagem em relação aos extremos
figure
imagesc(BW)
hold on

for i=1:l
    MamaLinha = Mama(i,:);
    valores = find(MamaLinha(1,:)>25 & MamaLinha(1,:)<38); %Filtrar os val-
% ores que correspondem a temperatura esperada para o corpo da pessoa, e
% não para o ambiente de fundo.

    MamaLinhaBW = BW(i,:);
    valoresBW = find(BW(i,:)==1);
    media(i,2) = i;
    media(i,1) = mean(valores); %Cálculo do ponto central da pessoa
    
    scatter(round(media(i,1)),round(media(i,2))); %Desenhar na imagem os 
% pontos de interesse.
end

%% Desconsiderar outliers referent aos pontos médios

% % Tentativa em se utilizar a moda e o desvio padrão para se retirar os
% outliers. Verificou-se que não funcionava por excluirem pontos de
% interesse.

% moda = mode(media(:,1));
% desvio = std(media(:,1));
% media(media(:,1)>moda-desvio & media(:,1)<moda+desvio,:)=[]

% Exclui-se aqui a região acima do seio visto que apresenta uma grande
% variação do centro devido a presença dos braços que raramente estão
% simétricos.

media(media(:,2)<250,:)=[];

% media(media(:,1)>moda-desvio & media(:,1)<moda+desvio,:)=[];

   [l c] = size(media);

% for i=1:l
% %     MamaLinha = Mama1(i,:);
% %     valores = find(MamaLinha(1,:)>25 & MamaLinha(1,:)<38);
% %     media(i,2) = i;
% %     media(i,1) = mean(valores);
%     scatter(round(media(i,1)),round(media(i,2)));
% end

%% Criar uma linha de simetria

mediaCorr = media(any(~isnan(media(:,1)),2),:); % Aqui se exclue os valores
% NaN que correspondem a valores não numéricos encontrados na média.

Coef = polyfit(mediaCorr(:,2),mediaCorr(:,1), 1);%Identifica-se qual a reta
% que mais se aproxima dos pontos de interesse.

Coef(:,1)=-Coef(:,1); %Correção da reta de modo a deixa-la na vertical.
x = 1:l; 
f = polyval(Coef,x)-l*Coef(:,1); %Deslocamento horizontal da reta.

plot(f,x,'-'); %Plotar a reta na imagem

%% Criar nova imagem rotacionada

B = imrotate(Mama,-atand(Coef(1))); %Rotação da imagem de acordo com a in-
% clinação da reta.
RotateBW = imrotate(BW,-atand(Coef(1))); % Rotação da imagem binaria de a-
%cordo com a inclinação da reta.

%   Plotar a nova imagem obtida.
    figure
    imagesc(B)
    hold on
    
% Obter o tamanho da nova imagem
[l c] = size(B);
for i=1:l
    MamaLinha = B(i,:);
    valores = find(MamaLinha(1,:)>25 & MamaLinha(1,:)<38); %Filtrar os val-
% ores que correspondem a temperatura esperada para o corpo da pessoa, e
% não para o ambiente de fundo.
    
    mediaN(i,2) = i;
    mediaN(i,1) = mean(valores); %Cálculo do ponto central da pessoa
    
    scatter(round(mediaN(i,1)),round(mediaN(i,2))); %Desenhar na imagem os 
% pontos de interesse.
end

%% Separar as imagens em duas simetricas e invertida

% Adotar a posição da reta para se dividir a imagem.
% OBS. Aqui evita-se que o valor seja negativo ou muito elevado.

centro=-l*Coef(:,1)/2+Coef(2); %Deslocamento da imagem para se realizar
% uma sobreposição adequada.

if centro<c/2 %Coef(2)
    colunas = centro;
else
    colunas = c-centro; 
end

% Criação e rotação das novas 2 imagens 
for i=1:l-1
    for j=1:colunas-1
    ImagemA(i,j) = B(i,round(centro+j));
    ImagemB(i,j) = B(i,round(centro-j));
    NovaImage(i,j) = abs(ImagemA(i,j)-ImagemB(i,j));
    end
end

figure
imagesc(NovaImage);

%% Deteção dos picos - Não operacional
% A detecção é importante para se saber o quão errado esta a sobreposição
% das imagens e tambem posteriormente para uma possivel detecção do cancer
% de mama.

% Não se adotou ainda o método de vargas devido a dificuldade de
% sobrepor-se adequadamente as imagens. Trabalhos futurão devem contemplar
% o método.

for i=1:size(NovaImage,2)-1;
Saida(i) = NovaImage(300,i)-NovaImage(300,i+1);
end
[XX YY] = findpeaks(abs(Saida),'NPeaks', 5,'MinPeakProminence',1 );

%% Criacão de uma lista por meio das bordas - Não operacional

% % A criação da lista é importante para se saber quais são os features
% para a detecção do pontos imutaveis na imagem 

% [x,y]=meshgrid(1:size(ImagemA,1), 1:size(ImagemA,2));
% ListaA=[x(:),y(:),ImagemABW(:)];
% ListaA(ListaA(:,3)==0,:)=[];
% ListaA(:,3)=[];
% 
% [x,y]=meshgrid(1:size(ImagemB,1), 1:size(ImagemB,2));
% ListaB=[x(:),y(:),ImagemBBW(:)];
% ListaB(ListaB(:,3)==0,:)=[];
% ListaB(:,3)=[];
%% Remover bordas da imagem - Incompleto
% Remover as linhas que possuem zero na imagem. Espera-se que as linhas
% removidas não possuam nenhuma informação relevante referente ao cancer.

% ImagemA1 = ImagemA;
% ImagemA(ImagemA(:,1)==0,:)=[];
% ImagemA(:,ImagemA(1,:)==0)=[];
% ImagemA(:,ImagemA(end,:)==0)=[];
% % 
% ImagemB(ImagemB(:,1)==0,:)=[];
% ImagemB(ImagemB(:,end)==0,:)=[];
% ImagemB(:,ImagemB(end,:)==0)=[];

%% Criar Lista por feature - Não operacional
% Espera-se utilizar o detector Harris para a identificação dos pontos de
% interesse.

cornersA = detectHarrisFeatures(ImagemA, 'MinQuality',0.00001);
ListaA = cornersA.Location;

cornersB = detectHarrisFeatures(ImagemB, 'MinQuality',0.00001)
ListaB = cornersB.Location;
% plot(corners.selectStrongest(10));

%% Criar listas com o mesmo tamanho
% Deixar as listas com o mesmo tamanho

if size(ListaA,1)-size(ListaB,1)<0
    diftam=size(ListaA,1)-size(ListaB,1);
    ListaBB = ListaB(abs(diftam)+1:end,:);
    ListaAA = ListaA;
else
    diftam=size(ListaA,1)-size(ListaB,1);
    ListaAA = ListaA(abs(diftam)+1:end,:);
    ListaBB = ListaB;
end

%% Calcular transformada adotada - Não operacional
%Identificar qual a transformada linear adequada para sobrepor
%adequadamente a imagem de modo qua a imagem seja esticada, rotacionada
%adequadamente para a sobreposição.

tform = fitgeotrans(ListaBB, ListaAA,'similarity');
% tform = estimateGeometricTransform(ListaBB,ListaAA,'projective');
ImagemBB = imwarp(ImagemB,tform);

% % Futuramente tentar-se-a utilizar as bordas para a detecção das
% % features.

% cornersBW = detectHarrisFeatures(BW, 'MinQuality',0.001);
% ListaB = cornersBW.Location;
figure
imagesc(ImagemA); hold on;
scatter(ListaAA(:,1),ListaAA(:,2));
figure
imagesc(ImagemB); hold on;
% plot(cornersBW.selectStrongest(100));
scatter(ListaBB(:,1),ListaBB(:,2));
%% Match Points - Não operacional

% % Apresentar na imagem os pontos de interesse e o quanto eles foram
% delocados.

% Espera-se utilizar o detector Harris para a identificação dos pontos de
% interesse.
points1 = detectHarrisFeatures(ImagemA);
points2 = detectHarrisFeatures(ImagemB);

[features1,valid_points1] = extractFeatures(ImagemA,points1);
[features2,valid_points2] = extractFeatures(ImagemB,points2);
% % O uso das Surf Features não se demonstrou funcional

% points1 = detectSURFFeatures(ImagemA,'MetricThreshold',10,'NumScaleLevels',6);
% points2 = detectSURFFeatures(ImagemB,'MetricThreshold',10,'NumOctaves',1);

indexPairs = matchFeatures(features1,features2,'MaxRatio',0.8);

matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

figure; showMatchedFeatures(ImagemA,ImagemB,matchedPoints1,matchedPoints2);
% legend('matched points 1','matched points 2');
