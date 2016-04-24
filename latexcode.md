# Double-Pendulum

\documentclass[12pt,a4paper,oneside]{book}
\usepackage[utf8]{inputenc}
\usepackage[hungarian]{babel}
\usepackage{amsmath,amsfonts,amssymb,graphicx}
\usepackage[justification=centering]{caption}
\usepackage{array}
\usepackage{subfig}
\usepackage{gensymb}
\usepackage{fancyhdr}
\usepackage{hyperref}			% hivatkozásokhoz
\hypersetup{colorlinks=false}	% hivatkozások színét lehet itt változtatni
\usepackage{indentfirst}		% az első bekezdés is behúzott legyen, ne csak a 2-tól
\linespread{1.3}				% másfeles sorköz
%\usepackage[none]{hyphenat}	% elválasztás kikapcsolása
\sloppy							% szóközök növelése szebb sorok érdekében
\usepackage[left=3.00cm, right=2.00cm, top=2.00cm, bottom=2.50cm]{geometry}

\begin{document}

\pagestyle{fancy}
\fancyhf{}

\chapter{Bevezetés}
%PÁRATLAN OLDALON!

\lhead{ \textit{1. Fejezet - Bevezetés}}
\rhead{ \textit{ \thepage}}
	
\chapter{Szakirodalmi áttekintés}
%PÁRATLAN OLDALON!

\lhead{ \textit{2. Fejezet - Szakirodalmi áttekintés}}
\rhead{ \textit{ \thepage}}
	
\chapter{A kettős inga mechanikai modellje}
%PÁRATLAN OLDALON!
	
\lhead{ \textit{3. Fejezet - A kettős inga mechanikai modellje}}
\rhead{ \textit{ \thepage}}

Ebben a fejezetben bemutatásra kerül a dolgozatban vizsgált kettős inga szerkezete, konstrukciós modellje. Bevezetjük a modellel kapcsolatos egyszerűsítéseket, feltételezéseket. Emellett megismerkedhetünk a rendszert leíró differenciálegyenletekkel, mozgásegyenletekkel, valamint azok meghatározásának menetével. A felírt, nemlineáris differenciálegyenletekből levezetésre kerülnek a lineáris mozgásegyenletek is, melyek később a rendszervizsgálat és a szabályozás alapjául fognak szolgálni. Végül néhány eseten keresztül demonstráljuk az adott kezdeti pozícióból elengedett, magára hagyott, szabadon lengő rendszer mozgását, viselkedését.

\begin{section}{A kettős inga szerkezeti modellje}

A \ref{szerkezetiabra} ábrán látható a kettős inga síkbeli szerkezeti ábrája. A kép egyensúlyi, nyugalmi helyzetében szemlélteti a szabályozatlan rendszert. A szerkezeti modell az $1$-es, illetve $2$-es jelöléssel ellátott homogén, prizmatikus, merev rudakból, továbbá a $3$-as jelölésű, egyenáramú motorból és a $4$-es alátámasztásból épül fel.
\begin{figure}[h]
\centering
\includegraphics[width=0.293\linewidth]{szerkezetiabra2}
\caption{A szabályozatlan kettős inga stabil, nyugalmi helyzetben}
\label{szerkezetiabra}
\end{figure}\par
A rudak, azaz inga-elemek anyaga alumínium, geometriájuk és méretük azonos, keresztmetszetük téglalap kialakítású.\par
A rendszer paramétereit tekintve legyen a belső, narancssárga elem tömege $m_{1}$, hossza pedig $l_{1}$. Hasonlóképpen legyen a külső, lila színnel jelzett tag tömege $m_{2}$, hossza $l_{2}$. Az itt bevezetett jelöléseket a \ref{Lagrangeabra}-es ábra is demonstrálja. A későbbiekben előforduló mennyiségek, paraméterek indexelése is ezen gondolatmenetet fogja követni.\par
Számszerűsítve az adatokat, az inga-elemek hossza $l_{1}=l_{2}=0.4\ [m]$, valamint tömege $m_{1}=m_{2}=\rho \cdot A_{1} \cdot l_{1}=\rho \cdot A_{2} \cdot l_{2}=0.108\ [kg]$. Ahol $A_{1}=A_{2}=20\ [mm] \times 5\ [mm]=10^{-4}\ [m^{2}]$ a keresztmetszet profilja, illetve $\rho=2700\ [kg/m^{3}]$ az alumínium szobahőmérsékleten vett sűrűsége [1!!!].\par
A kettős ingát alkotó rudak egy ideális, súrlódásmentes csuklón keresztül kapcsolódnak egymáshoz. A csapágyazás a motor hajtótengelyének és a belső inga-elemnek a csatlakozásánál szintén ideális, súrlódásmentes. A közegellenállás által a szerkezetre ható erő elhanyagolható, így a modellezés során ezt a hatást is figyelmen kívül hagyjuk. Mindezekből adódóan megállapítható, hogy a szabályozatlan inga nem veszít energiát mozgása, lengése során, hiszen a rendszer nem tartalmaz disszipatív mennyiségeket. Az elemek mozgási és helyzeti energiái különböző módon, veszteségmentesen alakulnak át egymásba.\par
Továbbá, mivel a dolgozatnak nem célja a szerkezet elektronikai egységeinek modellezése vagy vizsgálata, a rendszert leíró differenciálegyenletek nem tartalmazzák a DC motor dinamikáját, működését jellemző egyenleteket. A mechanikai modell meghatározásakor a motor által a rendszerre kifejtett nyomatékot pusztán egy nyomatékgerjesztéssel helyettesítjük.
\end{section}

\begin{section}{A rendszer mozgásegyenleteinek meghatározása}

\begin{subsection}{A másodfajú-Lagrange egyenlet}

Egy adott mechanikai modell mozgásegyenleteinek felírására igen gyakran alkalmazott módszer az úgynevezett Newton-Euler módszer. Ekkor a szükséges számú szabadtest ábra felrajzolása után, a dinamika alaptételét alkalmazva juthatunk el a keresett differenciálegyenletekhez. Azonban összetettebb, több szabadságfokú mechanikai rendszerek esetén ez az eljárás hosszadalmas lehet. Ilyenkor célravezetőbb úgymond analitikus úton, általában energia- vagy teljesítmény jellegű mennyiségek meghatározása és differenciálása révén felírni a mozgásegyenleteket. Az egyik legfontosabb és legelterjedtebb ilyen analitikus módszer a teljesítmény tétel továbbfejlesztésének tekinthető másodfajú Lagrange-egyenlet alkalmazása [2!!!].\par
A Lagrange-módszer használata során elsőként intuitív módon meg kell választanunk az egymástól független általános koordinátákat a rendszer szabadsági fokának megfelelően. Majd ezután a kiválasztott általános koordináták segítségével fel kell írnunk a rendszer kinetikus energiáját, valamint az aktív erőket és a nem ideális kényszererőket tartalmazó, úgynevezett általános erő vektort. Mindezek felhasználásával a másodfajú Lagrange-egyenlet az alábbi alakot ölti [2!!!]:
\begin{equation}
\frac{\textup{d}}{\textup{d} t} \frac{\partial E_{K}}{\partial \dot{q_{j}}} - \frac{\partial E_{K}}{\partial q_{j}} = Q_{j}, \qquad  j=1,...n.
\label{alapLagrange}
\end{equation}
Ahol $n$ a rendszer szabadsági fokainak száma, $q_{j}$ a $\mathbf{q}=(q_{1},q_{2},\dots q_{n})^{T}$ általános koordináta vektor $j$-edik eleme, $E_{K}$ a rendszer kinetikus energiája, illetve $Q_{j}=\sum\limits_{i=1}^N \mathbf{F_{i}} \frac{\partial \mathbf{r_{i}}}{\partial q_{j}}$ az általános erő vektorának $j$-edik komponense. Itt $\mathbf{F_{i}}$ az $i$-edik tömegpontra ható aktív erőket és a nem ideális kényszererőket jelöli, $\mathbf{r_{i}}$ az $i$-edik tömegpont helyvektora az általános koordinátákkal kifejezve, $N$ a tömegpontok száma [2!!!].\par
A \ref{alapLagrange} összefüggés megadja egy $n$ szabadsági fokú mechanikai rendszer mozgásegyenleteit. Ugyanakkor, ha az általános erő vektorának komponenseit felbontjuk potenciálos erőkre (pl. nehézségi erő, rugóerő), csillapító erőkre (pl. súrlódás, viszkózus csillapítás) és nem potenciálos erőkre, akkor a másodfajú Lagrange-egyenletnek a következő, gyakorlatban jobban használható alakját kapjuk [2!!!]:
\begin{equation}
\frac{\textup{d}}{\textup{d} t} \frac{\partial E_{K}}{\partial \dot{q_{j}}} - \frac{\partial E_{K}}{\partial q_{j}} + \frac{\partial \mathcal{D}}{\partial \dot{q_{j}}} + \frac{\partial U}{\partial q_{j}} = Q_{j}^*, \qquad  j=1,...n.
\label{Lagrange}
\end{equation}
Ahol $U$ a rendszer potenciális energiája, $\mathcal{D}$ a disszipatív potenciál, továbbá $Q_{j}^*$ az általános erő $j$-edik komponensének az a része, ami nem fejezhető ki az egyenlet bal oldalán szereplő mennyiségek segítségével [2!!!].
\end{subsection}

\begin{subsection}{A kettős inga nemlineáris mozgásegyenletei}
Célszerű a Lagrange-módszert alkalmazni a kettős inga mozgásegyenleteinek meghatározásához is. Ehhez először tekintsük az alábbi, \ref{Lagrangeabra} ábrát.
\begin{figure}[h]
\centering
\includegraphics[scale=0.1269]{szerkezetiabra4}
\caption{A szabályozatlan kettős inga felső, instabil egyensúlyi helyzetéből kitérítve}
\label{Lagrangeabra}
\end{figure}\par
A képen látható $(x,y)$ koordináta-rendszer mutatja a pozitív elmozdulási irányokat. A szögelfordulások esetében az óramutató járásával ellentétes irányt tekintjük pozitívnak.
A mechanikai rendszerünk két szabadsági fokú ($n=2$), így két egymástól független általános koordinátát szükséges választanunk. Ahogy \ref{Lagrangeabra} ábra is szemlélteti, kézenfekvő az inga-elemek felső, függőleges helyzetétől mért szögkitéréseket általános koordinátáknak választani, azaz: $q_{1}=\Phi_{1}$ és $q_{2}=\Phi_{2}$. Vagyis az általános koordináták vektora:
\begin{equation}
\mathbf{q}=\begin{bmatrix}
\Phi_{1} \\
\Phi_{2}
\end{bmatrix}
\label{altkoord}
\end{equation}\par
$\mathbf{q}$ vektor ismeretében nincs más dolgunk, minthogy a rendszerre vonatkoztatva felírjuk a \ref{Lagrange} egyenletben szereplő mennyiségeket.
Mivel az inga a nehézségi erőtérben helyezkedik el, kinetikus energián kívül potenciális energiával is rendelkezik. A disszipatív potenciál azonban zérus, hiszen nincsen olyan hatás, ami a rendszer energiáját csökkenteni igyekszik.\\ [0.2 cm]
A belső inga-elem kinetikus energiája: \vspace{1 mm} 
\begin{equation} \vspace{1 mm} 
E_{K,1}=\frac{1}{2} m_{1} v_{S1}^{2}+\frac{1}{2} \mathit{\Theta}_{S1} \dot{\mathit{\Phi}}_{1}^{2}=\frac{1}{6} m_{1} l_{1}^{2} \dot{\mathit{\Phi}}_{1}^{2}
\label{ek1}
\end{equation}
Ahol $\underline{v}_{S1}=\underline{v}_{0} + \dot{\mathit{\underline{\Phi}}}_{1} \times \underline{r}_{0S1}=(- \frac{l_{1}}{2} \cos(\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1},\ - \frac{l_{1}}{2} \sin(\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1})^{T}$ a belső inga súlypontjának sebességvektora. Itt $\underline{r}_{0S1}$ az origóból az $S1$ súlypontba mutató helyvektor, illetve $\underline{v}_{0}$ a belső csukló sebességvektora. Mivel a belső csukló - rögzítettségéből adódóan - sebességpólus, ezért  $\underline{v}_{0}=\underline{0}$. Továbbá $\mathit{\Theta}_{S1} \cong \frac{1}{12} m_{1} l_{1}^2$ az inga-elemnek a saját súlypontjára számított tehetetlenségi nyomatéka, valamint $\dot{\mathit{\underline{\Phi}}}_{1}=\underline{\omega}_{1}$ a belső inga szögsebesség-vektora. \\ [0.2 cm]
A külső inga-elem kinetikus energiája: \vspace{1 mm}
\begin{equation} \vspace{1 mm}
E_{K,2}=\frac{1}{2} m_{2} v_{S2}^{2}+\frac{1}{2} \mathit{\Theta}_{S2} \dot{\mathit{\Phi}}_{2}^{2}=\frac{1}{2} m_{2} l_{1}^{2} \dot{\mathit{\Phi}}_{1}^{2} + \frac{1}{6} m_{2} l_{2}^{2} \dot{\mathit{\Phi}}_{2}^{2} + \frac{1}{2} m_{2} l_{1} l_{2} \cos(\mathit{\Phi}_{2}-\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1} \dot{\mathit{\Phi}}_{2}
\label{ek2}
\end{equation}
Ahol $\underline{v}_{S2}=\underline{v}_{1} + \dot{\mathit{\underline{\Phi}}}_{2} \times \underline{r}_{1S2}=(-l_{1} \cos(\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1} - \frac{l_{2}}{2} \cos(\mathit{\Phi}_{2}) \dot{\mathit{\Phi}}_{2},\ - l_{1} \sin(\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1} - \frac{l_{2}}{2} \sin(\mathit{\Phi}_{2}) \dot{\mathit{\Phi}}_{2})^{T}$ a külső inga súlypontjának sebességvektora. Itt $\underline{r}_{1S2}$ a külső csuklóból az $S2$ súlypontba mutató helyvektor, illetve $\underline{v}_{1}=2 \underline{v}_{S1}$ a külső csukló sebességvektora. Továbbá $\mathit{\Theta}_{S2} \cong \frac{1}{12} m_{2} l_{2}^2$ a külső inga-elemnek a saját súlypontjára számított tehetetlenségi nyomatéka, valamint $\dot{\mathit{\underline{\Phi}}}_{2}=\underline{\omega}_{2}$ a külső inga szögsebesség-vektora.\\ [0.2 cm]
Így a rendszerünk teljes kinetikus energiája: \vspace{1 mm}
\begin{equation} \vspace{1 mm}
E_{K}=E_{K,1}+E_{K,2}=\left(\frac{1}{6} m_{1} + \frac{1}{2} m_{2} \right) l_{1}^{2} \dot{\mathit{\Phi}}_{1}^{2} +  \frac{1}{6} m_{2} l_{2}^{2} \dot{\mathit{\Phi}}_{2}^{2} + \frac{1}{2} m_{2} l_{1} l_{2} \cos(\mathit{\Phi}_{2}-\mathit{\Phi}_{1}) \dot{\mathit{\Phi}}_{1} \dot{\mathit{\Phi}}_{2}
\label{ekfull}
\end{equation}
A belső inga-elem potenciális energiája:
\begin{equation}
U_{1}=\frac{1}{2} m_1 g l_1 \cos(\mathit{\Phi}_{1})
\label{ep1}
\end{equation}
Ahol $g=9.81\ [m/s^2]$ a nehézségi gyorsulás nagysága. \\ [0.2 cm]
A külső inga-elem potenciális energiája:
\begin{equation} 
U_{2}=m_2 g \left(l_1 \cos(\mathit{\Phi}_{1}) + \frac{1}{2} l_2 \cos(\mathit{\Phi}_{2}) \right)
\label{ep2}
\end{equation}
Tehát a rendszer teljes potenciális energiája:
\begin{equation}
U=U_{1} + U_{2}=\frac{1}{2} m_1 g l_1 \cos(\mathit{\Phi}_{1}) + m_2 g \left(l_1 \cos(\mathit{\Phi}_{1}) + \frac{1}{2} l_2 \cos(\mathit{\Phi}_{2}) \right)
\label{epfull}
\end{equation}
A disszipatív potenciál, ahogy már említésre került:
\begin{equation}
\mathcal{D}=0
\label{disp}
\end{equation}
Végül definiáljuk a másodfajú Lagrange-egyenlet jobb oldalán szereplő $\mathbf{Q}^*$ vektort. Az egyetlen hatás, amit még nem vettünk figyelembe a fent meghatározott mennyiségeken keresztül, az a nyomatékgerjesztés, azaz a motor által a belső inga-elemre kifejtett nyomaték. Így a $\mathbf{Q}^*$ vektor az alábbi alakot ölti:
\begin{equation}
\mathbf{Q}^*=\begin{bmatrix}
	M \\
	0
\end{bmatrix}
\label{Q}
\end{equation}

A kettős inga mozgásegyenleteinek felírásához már csak származtatnunk kell a fenti mennyiségeknek a \ref{Lagrange} egyenlet szerinti deriváltjait. A deriváltak előállítása után a rendszert leíró két nemlineáris differenciálegyenlet:

\begin{equation}
\label{eq1}
\begin{gathered}
\left(\frac{1}{3} m_1 + m_2 \right)l_1^2 \ddot{\mathit{\Phi}}_{1} + \frac{1}{2} m_2 l_1 l_2 \cos(\mathit{\Phi}_2-\mathit{\Phi}_1) \ddot{\mathit{\Phi}}_{2} - \frac{1}{2} m_2 l_1 l_2 \sin(\mathit{\Phi}_2-\mathit{\Phi}_1) \dot{\mathit{\Phi}}_{2}^2 \\
- \left(\frac{1}{2} m_1 + m_2 \right) g l_1 \sin(\mathit{\Phi}_{1})=M
\end{gathered}
\end{equation}

\begin{equation} \vspace{5 mm}
\label{eq2}
\frac{1}{2} m_2 l_1 l_2 \cos(\mathit{\Phi}_2-\mathit{\Phi}_1) \ddot{\mathit{\Phi}}_{1} + \frac{1}{3} m_2 l_2^2 \ddot{\mathit{\Phi}}_{2} + \frac{1}{2} m_2 l_1 l_2 \sin(\mathit{\Phi}_2-\mathit{\Phi}_1) \dot{\mathit{\Phi}}_{1}^2 - \frac{1}{2} m_2 g l_2 \sin(\mathit{\Phi}_{2})=0
\end{equation}

\end{subsection}

\begin{subsection}{A kettős inga lineáris mozgásegyenletei}
A kettős ingát a $\mathit{\Phi}_{1}=\mathit{\Phi}_{2}=0\ [\degree]$, felső, függőleges helyzetében szeretnénk stabilizálni egy megfelelően megtervezett lineáris szabályozás segítségével. Ezen lineáris szabályozás megvalósításához először is meg kell határoznunk a rendszert leíró lineáris differenciálegyenleteket. Feltételezve, hogy az inga felső, nulla pont körüli kitérései megfelelően kicsik lesznek a szabályozás során, linearizálhatjuk a \ref{eq1} és \ref{eq2} mozgásegyenleteket.\par
Linearizálás esetén az adott nemlineáris egyenletet, vagy az egyenlet nemlineáris tagjait az elsőrendű Taylor-polinommal közelítjük. Tehát az elsőfokú közelítést úgy kapjuk, hogy a kifejezések Taylor-sorából elhagyjuk a másod- vagy magasabb fokú tagokat [2] [3!!].\par
Megvizsgálva  a \ref{eq1} és \ref{eq2} differenciálegyenleteket, láthatjuk, hogy a nemlineáris tagok vagy a szögelfordulások egyváltozós, trigonometrikus függvényei, vagy a szögelfordulások mellett a szögsebességeket, illetve szöggyorsulásokat is tartalmazó háromváltozós függvények. Ebből adódóan az egy- és a háromváltozós függvényekre érvényes elsőrendű Taylor-polinomok szolgáltatják a keresett lineáris kifejezéseket. \\[0.2 cm]
Az egyváltozós $f(x)$ függvény $x=0$ pontbeli elsőrendű Taylor-polinomja [3]:
\begin{equation}
f(x) \cong f(0) + f^{'}(0) \cdot x
\label{tp1}
\end{equation}
A háromváltozós $f(x,y,z)$ függvény $(x,y,z)=(0,0,0)$ pontbeli elsőrendű Taylor-polinomja [3]:
\begin{equation}\vspace{2 mm}
f(x,y,z) \cong f(0,0,0) + \left. \left(\frac{\partial f(x,y,z)}{\partial x} \cdot x + \frac{\partial f(x,y,z)}{\partial y} \cdot y + \frac{\partial f(x,y,z)}{\partial z} \cdot z \right)\right \vert_{(0,0,0)}
\label{tp3}
\end{equation}
A \ref{tp1} és \ref{tp3} összefüggéseket felhasználva az egyes nemlineáris tagok linearizált alakja:
\begin{itemize}
\item $\sin(\mathit{\Phi}_{1}) \cong \mathit{\Phi}_{1} \qquad \qquad \qquad \sin(\mathit{\Phi}_2-\mathit{\Phi}_1) \dot{\mathit{\Phi}}_{2}^2 \cong 0 \qquad \qquad \qquad \cos(\mathit{\Phi}_2-\mathit{\Phi}_1) \ddot{\mathit{\Phi}}_{2} \cong \ddot{\mathit{\Phi}}_{2}$
\item $\sin(\mathit{\Phi}_{2}) \cong \mathit{\Phi}_{2} \qquad \qquad \qquad \sin(\mathit{\Phi}_2-\mathit{\Phi}_1) \dot{\mathit{\Phi}}_{1}^2 \cong 0 \qquad \qquad \qquad \cos(\mathit{\Phi}_2-\mathit{\Phi}_1) \ddot{\mathit{\Phi}}_{1} \cong  \ddot{\mathit{\Phi}}_{1}$
\end{itemize}\par \vspace{2 mm}
A kettős inga linearizált mozgásegyenleteit megkapjuk, ha a \ref{eq1} és \ref{eq2} egyenleteken belül az eredeti, nemlineáris tagok helyére behelyettesítjük a fent meghatározott linearizált kifejezéseket:

\begin{equation}
\label{leq1}
\left(\frac{1}{3} m_1 + m_2 \right)l_1^2 \ddot{\mathit{\Phi}}_{1} + \frac{1}{2} m_2 l_1 l_2 \ddot{\mathit{\Phi}}_{2} - \left(\frac{1}{2} m_1 + m_2 \right) g l_1 \mathit{\Phi}_{1}=M
\end{equation}

\begin{equation} \vspace{5 mm}
\label{leq2}
\frac{1}{2} m_2 l_1 l_2 \ddot{\mathit{\Phi}}_{1} + \frac{1}{3} m_2 l_2^2 \ddot{\mathit{\Phi}}_{2} - \frac{1}{2} m_2 g l_2 \mathit{\Phi}_{2}=0
\end{equation}\par

A kapott lineáris differenciálegyenletek kis kitérések esetén megfelelő pontossággal közelítik az eredeti mozgásegyenleteket, így felhasználhatjuk őket a rendszer szabályozása során.
\end{subsection}
\end{section}

\begin{section}{A kettős inga néhány mozgástörvénye}

Adott kezdeti feltételek mellett megoldva a mechanikai rendszerünket leíró mozgásegyenleteket, megkapjuk a kettős inga egy adott mozgástörvényét. Abban az esetben, ha a szabályozatlan ingának az alsó, egyensúlyi helyzete körüli kis kitéréseit akarjuk vizsgálni, akkor az ezen egyensúlyi helyzet körül linearizált egyenletek is kielégítő megoldást adhatnak a rendszer válaszát illetően. Azonban, ha a nagyobb kezdeti szögsebességekkel rendelkező vagy nagyobb szögben kitérített, akár a felső pozíció környezetéből elengedett, magára hagyott inga mozgását szeretnénk analizálni, akkor csak az eredeti, nemlineáris differenciálegyenletek szolgáltatnak megfelelő eredményt.\par
A \ref{eq1} és \ref{eq2} nemlineáris egyenleteknek ugyanakkor nincs zárt alakú megoldása, így numerikus úton, Wolfram Mathematica segítségével fogjuk meghatározni az adott kezdeti pozícióból elengedett és adott kezdeti szögsebességekkel rendelkező kettős inga mozgástörvényét.\\
Az alábbi ábrák két különböző kezdeti feltétel esetén szemléltetik a $\mathit{\Phi}_{1}$ és $\mathit{\Phi}_{2}$ szögelfordulás függvényeket.

\begin{figure}[h]
	\centering
	\subfloat[A belső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi11}}} 
	\,
	\subfloat[A külső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi21}}}
	\caption{A mozgástörvény $\mathit{\Phi}_{1}(0)=-140\ [\degree]$, $\mathit{\Phi}_{2}(0)=-150\ [\degree]$ és ${\omega}_{1}(0)={\omega}_{2}(0)=0\ [rad/s]$ kezdeti feltételek esetén}
	\label{mozgtorv1}
\end{figure} \vspace{-6 mm}
\begin{figure}[h]
	\centering
	\subfloat[A belső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi12}}} 
	\,
	\subfloat[A külső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi22}}}
	\caption{A mozgástörvény $\mathit{\Phi}_{1}(0)=-115\ [\degree]$, $\mathit{\Phi}_{2}(0)=-75\ [\degree]$ és ${\omega}_{1}(0)={\omega}_{2}(0)=0\ [rad/s]$ kezdeti feltételek esetén}
	\label{mozgtorv2}
\end{figure}

A magára hagyott kettős inga az egyensúlyi pontja ($\mathit{\Phi}_{1}=\mathit{\Phi}_{2}=-180\ [\degree]$) körül csillapítatlan rezgőmozgást végez. Már ezen két eseten keresztül is megállapíthatjuk, hogy a mozgás képe meglehetősen összetett és szabálytalan.\par
A rendszerünk érzékeny a kezdeti feltételekre. Minél nagyobb a kezdeti mozgási-, illetve helyzeti energia, azaz minél nagyobb a rendszer összenergiája, annál változatosabb, komplexebb a mozgástörvény. Különösen igaz ez akkor, amikor az egyes inga-elemek többször is átfordulnak a csuklópontjuk körül. Egy ilyen esetet szemléltet a \ref{mozgtorv3} ábra az adott felső pozícióból elengedett inga feltűnően szabálytalan mozgástörvényének segítségével. \par

\begin{figure}[h]
	\centering
	\subfloat[A belső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi13}}} 
	\,
	\subfloat[A külső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.49\linewidth]{szabalyozatlanfi23}}}
	\caption{A mozgástörvény $\mathit{\Phi}_{1}(0)=-20\ [\degree]$, $\mathit{\Phi}_{2}(0)=-10\ [\degree]$ és ${\omega}_{1}(0)={\omega}_{2}(0)=0\ [rad/s]$ kezdeti feltételek esetén}
	\label{mozgtorv3}
\end{figure}

Lent a \ref{mozgtorv1}, valamint a \ref{mozgtorv2} esetekhez tartozó trajektóriák $(\mathit{\Phi}_{1}$-$\mathit{\Phi}_{2})$ síkra való vetületeit figyelhetjük meg.
\begin{figure}[h]
	\centering
	\subfloat[\ref{mozgtorv1} esethez tartozó trajektória vetülete]{{\includegraphics[width=0.48\linewidth]{szabalyozatlantraj1}}} 
	\qquad \qquad
	\subfloat[\ref{mozgtorv2} esethez tartozó trajektória vetülete]{{\includegraphics[width=0.355\linewidth]{szabalyozatlantraj2}}}
	\caption{A trajektóriák $(\mathit{\Phi}_{1}$-$\mathit{\Phi}_{2})$ síkokra való vetületei}
	\label{traj}
\end{figure}

A trajektóriák vetületei is komplex, de ugyanakkor rendezett struktúrájú mozgásra utalnak. Ez a tény, valamint a kezdeti feltételekre való érzékenység, a hosszú távú aperiodikusság, a mozgás időbeli szabálytalansága mind azt mutatják, hogy a csillapítatlan kettős inga kaotikus jellegű rezgőmozgást végez.

\end{section}

\chapter{A kettős inga fellendítése}
%PÁRATLAN OLDALON!

\lhead{ \textit{4. Fejezet - A kettős inga fellendítése}}
\rhead{ \textit{ \thepage}}

Ezen fejezeten belül említést teszünk néhány szabályozási eljárásról, amely alkalmas a kettős inga felső pozícióba juttatására. Ezt követően ismertetjük a fellendítéshez alkalmazott módszert, majd demonstráljuk az így kapott eredményeket.

\begin{section}{A kettős ingák fellendítésének néhány lehetséges módszere}

A kettős ingák fellendítésének - csakúgy, mint az instabil helyzetben történő szabályozásuknak - nincs olyan módszere, mely egyértelműen a legoptimálisabb, legmegfelelőbb volna. Mindegyik eljárásnak megvan a maga előnye vagy éppen hátránya, nehézsége. A számunkra legideálisabb megoldás más és más lehet a rendszer felépítésétől, paramétereitől, valamint attól függően, hogy mit tekintünk fontosnak a fellendítés során (minimális nagyságú beavatkozó hatás, minél gyorsabb fellendülés, minél kisebb helyzeti energia a felső pozícióban stb.).\par
A továbbiakban a teljesség igénye nélkül megemlítünk pár módszert, mely lehetővé teszi, hogy a kettős ingát az alsó, nyugalmi helyzetéből a felső, instabil állapot környezetébe juttassuk. Az egyes szabályozásokat nem ismertetjük részletesen, csupán azok típusáról, valamint a fellendítés elvéről ejtünk szót.\par
A [4!!!]-es irodalmi forrás egy kiskocsin lévő kettős inga fellendítését s egyben szabályozását írja le. A szabályozás három lépésből áll. Elsőként a belső ingát lendítjük fel az inverz ingák témakörében már jól ismert energiaszabályozás útján. Ezt követően, a második lépés során egy hibrid szabályozó ugyancsak energiaszabályozással fellendíti a külső inga-elemet úgy, hogy eközben egy nemlineáris szabályozás, úgynevezett csúszómód szabályozás (Sliding Mode Control) stabilizálja a belső tagot. Végül a harmadik, utolsó lépésben, csúszómód szabályozást alkalmazva stabilizálásra kerül a teljes inga a felső, instabil egyensúlyi helyzetében.\par
A belső csuklópontjában nyomatékkal gerjesztett kettős inga fellendítésére láthatunk egy lehetséges megoldást az [5!!!]-ös szakirodalmi publikációban. A rendszer egy bang-bang típusú, visszacsatolásos szabályozás hatására jut a felső pozíció kis környezetébe. A bang-bang szabályozásból adódóan a gerjesztő nyomaték értéke két állapot ($\pm M_{max}$) között váltakozik a visszacsatolt szögelfordulásoktól és szögsebességektől függően. A megfelelően előállított szabályozójel egyre jobban kitéríti az ingát az egyensúlyi helyzetéből, a maximális kitérés egyre nő, mígnem a rendszer elér a felső pozíció környezetébe.\\
A cikkben részletesen leírt módszer különösen alkalmas nagy tömegű, nagy tehetetlenségű ingák kis szabályozó nyomatékkal történő fellendítésére. Ugyanakkor az inga-elemek közti szögkülönbség viszonylag nagy lehet a felső pozíciónál, így egy lineáris szabályozó nem feltétlenül képes stabilizálni a rendszert.\par
A [6]-os forrás szintén egy nyomatékkal gerjesztett kettős inga fellendítésével kapcsolatban ismertet egy visszacsatolásos szabályozási eljárást. Itt azonban a rendszer úgy jut a felső állapotába, hogy közben az inga-elemek közötti szögkülönbséget szabályozzuk. Ebből adódóan elérhetjük, hogy a felső, instabil pozíció környezetben a szögkülönbség elegendően kicsi legyen ahhoz, hogy lineáris szabályozással stabilizálni tudjuk a rendszert. Ehhez azonban az is szükséges, hogy a felső helyzetben az inga kinetikus energiája kellően kicsi legyen.\par
Előszeretettel alkalmazzák a fent említett energia-, illetve bang-bang típusú szabályozásokat különféle ingák fellendítésére, de ezeken kívül számos más megoldással is találkozhatunk.
\end{section}
	
\begin{section}{A kettős inga fellendítésére alkalmazott módszer}

Az [5]-ös publikációban leírt módszer segítségével akár a mi rendszerünk is a felső, instabil egyensúlyi helyzet kis környezetébe juttatható. Ugyanakkor fontos megjegyezni, hogy a szerzők által vizsgált esetben a cél feltehetőleg az volt, hogy egy nagyobb tömegű, nagyobb tehetetlenségű lengőrendszert viszonylag kis maximális nyomatékú motor segítségével lendítsenek fel. Esetünkben az első számú szempont azonban az volna, hogy a kettős ingát olyan állapotában juttassuk a felső pozíció környezetébe, hogy azt egy lineáris szabályozó képes legyen egyensúlyban tartani. Ehhez egyrészt szükséges, hogy amikor a rendszer felér a felső, függőleges helyzet környezetébe, akkor az egyes inga-elemek közti szögkülönbség a lehető legkisebb legyen, másrészt, hogy a kinetikus energia is minimális legyen. Ha mindezek nem teljesülnek, akkor a szabályozás során a belső ingának túl nagy szögkitéréseket kell ahhoz tennie, hogy a rendszert stabilizálni tudjuk, s így a lineáris szabályozás érvényét veszti.\par 
Ahogy már korábban is említettük, az [5]-ben bemutatott módszer sok esetben nem teszi lehetővé a lineáris szabályozó alkalmazását. Ebből adódóan célszerű lenne egy másik, olyan fellendítési elvet találni, amely amellett, hogy ezt lehetővé teszi, gyorsan feljuttatja az ingát az instabil egyensúlyi helyzet környezetébe.\par
A mechanikai rendszerünk tehetetlensége meglehetősen kicsi, így egy kisebb nyomatékú motor is könnyen mozgásba hozza a szerkezetet. Megfelelő nagyságú nyomaték képes ellenhatni nemcsak az inga tehetetlenségének, de a gravitációs erőnek is, így képes körbeforgatni az ingát a belső csuklópont körül. Ebből adódik az az elgondolás, hogyha a belső csuklópontra csupán egy rövid ideig egy meghatározott, állandó nagyságú nyomatékot kapcsolunk, akkor a kettős inga fellendülhet. A belső csuklópontban ható konstans nyomaték elkezdi felfelé forgatni a belső ingát. A külső tag tehetetlenségének következtében kezdetben ugyan lemarad, de a belső inga-elem úgymond "magával rántja", sebessége, szögsebessége egyre nő, így utoléri a belső tagot, szögkülönbségük minimális lesz, tengelyeik szinte egybeesnek. A nyomaték kikapcsolása utána a meglévő kinetikus energia egy része fokozatosan átalakul helyzeti energiává, ahogy a két inga-elem együtt halad a felső, függőleges állapot felé. Adódhat olyan eset, hogy a mozgási energia jelentős hányada helyzetivé alakul a felső pozíció környezetében, így a kettős inga úgy éri el az instabil egyensúlyi helyzetet, hogy az egyes tagok sebessége, szögsebessége minimális, továbbá a szögkülönbségük is optimális.\par
A felvázolt módszer meglehetősen egyszerű, ugyanakkor esetünkben különösen kézenfekvő. Lehetővé teszi a kettős inga gyors fellendítését, továbbá azt, hogy a felső, instabil pozícióban egy lineáris szabályozó segítségével stabilizálni tudjuk a rendszert.\par
A módszer sikerének kulcsa, hogy találjunk egy, a rendszerünk paramétereinek megfelelő olyan nyomatékot, melyet, ha alkalmas időpillanatban lekapcsolunk, akkor az inga a számunkra ideális módon, úgy jut a felső pozíció környezetébe, hogy ott a szögkülönbség és a szögsebességek értéke is minimális. Ehhez legcélszerűbb, ha numerikus szimulációkat végzünk.
\end{section}

\begin{section}{A fellendítés szimulációjának eredményei}

A fellendítés szimulációja során vizsgáljuk, hogy mekkora nagyságú és mennyi ideig ható nyomaték képes a kettős ingát a stabil, nyugalmi helyzetéből a felső, instabil egyensúlyi helyzet környezetébe juttatni a már említett, ideális módon.\par 
A nyomatékot nem egy adott idő eltelte után kapcsoljuk le, hanem a belső inga egy megadott szögkitérése után. Abban az esetben, ha a belső inga-elem eléri a felső, függőleges pozíciót ($\mathit{\Phi}_1=0\ [\degree]$), akkor úgy tekintjük, hogy a kettős inga valamilyen állapotában fellendült. Ebben az állapotban, mikor $\mathit{\Phi}_1=0\ [\degree]$, vizsgáljuk az inga-elemek közti szögkülönbséget, a szögsebességeket, valamint a fellendítés idejét. Ezen adatokat mérlegelve választjuk majd ki, hogy mekkora legyen a nyomaték nagysága, illetve, hogy a belső inga mekkora szögkitérése utána kapcsoljuk le azt.\par
Fokozatosan növelve a gerjesztő nyomatékot, megtalálhatjuk azt a határértéket, amivel már képesek vagyunk valamilyen módon fellendíteni a rendszert. Ezt követően egy megadott tartományon belül, adott lépésközökkel változtatjuk a nyomaték nagyságát úgy, hogy egy-egy nyomatékértékhez tartozóan a kikapcsolási szög értékét is adott határok között, megadott lépésközzel változtatjuk. Miután találtunk az elvárásainkhoz közelítő értékeket, érdemes az adott nyomaték körül egy szűkebb tartományt venni, illetve finomítani a lépésközt. Hasonlóan a kikapcsolási szög értékét is vizsgálhatjuk szűkebb tartományon, kisebb lépésközökkel. Ezt az elvet követve olyan paraméterekhez juthatunk, melyek már teljes egészében kielégítik a fellendítéssel kapcsolatos követelményeinket.\par
Az alábbi két táblázat olyan eseteket szemléltet, ahol már egy konkrét nyomatékérték mellett, egy meglehetősen szűk tartományon belül, kis lépésközzel keressük az ideális kikapcsolási szöget. \\

\begin{table}[h]
	\begin{center}
		\begin{tabular}{ | >{\centering}m{2cm} | >{\centering}m{2.4cm}| >{\centering}m{2.4cm} | >{\centering}m{2cm} | >{\centering}m{2cm}| >{\centering}m{2cm} | }
			\hline
			\textit{M} [Nm] & \textit{Kikapcsolási szög} $[\degree]$ & \textit{Fellendítési idő} [s] & $\mathit{\Phi}_2-\mathit{\Phi}_1 [\degree]$ & $\mathit{\omega_1}$ [rad/s] & $\mathit{\omega_2}$ [rad/s]  \tabularnewline \hline
			\textbf{0.738} & \textbf{-48.00} & \textbf{1.23} & \textbf{-1.4904} & \textbf{0.4943} &  \textbf{0.4225}  \tabularnewline
			0.738 & -47.99 & 1.15 & -7.6372 & 0.8640 & -0.2268  \tabularnewline
			0.738 & -47.98 & 1.12 & -10.2990 & 1.0307 & -0.4739  \tabularnewline
			0.738 & -47.97 & 1.10 & -12.0952 & 1.1442 & -0.6270  \tabularnewline
			0.738 & -47.96 & 1.09 & -13.4735 & 1.2314 & -0.7367  \tabularnewline
			0.738 & -47.95 & 1.08 & -14.6012 & 1.3027 & -0.8211  \tabularnewline
			0.738 & -47.94 & 1.07 & -15.5604 & 1.3634 & -0.8890  \tabularnewline \hline
		\end{tabular}
	\end{center}
	\caption{$0.738$ [Nm] nagyságú nyomatékhoz tartozó adatok}
	\label{tablazat1}
\end{table}

\begin{table}[h]
	\begin{center}
		\begin{tabular}{ | >{\centering}m{2cm} | >{\centering}m{2.4cm}| >{\centering}m{2.4cm} | >{\centering}m{2cm} | >{\centering}m{2cm}| >{\centering}m{2cm} | }
			\hline
			\textit{M} [Nm] & \textit{Kikapcsolási szög} $[\degree]$ & \textit{Fellendítési idő} [s] & $\mathit{\Phi}_2-\mathit{\Phi}_1 [\degree]$ & $\mathit{\omega_1}$ [rad/s] & $\mathit{\omega_2}$ [rad/s] \tabularnewline \hline
			0.739 & -47.88 & 1.19 &  1.6587  & 0.4529 & 0.9841 \tabularnewline
			0.739 & -47.87 & 1.13 & -4.7700  & 0.8186 & 0.2634 \tabularnewline
			0.739 & -47.86 & 1.11 & -7.6579  & 0.9915 & -0.0313 \tabularnewline
			0.739 & -47.85 & 1.09 & -9.6129  & 1.1101 & -0.2174 \tabularnewline
			0.739 & -47.84 & 1.08 & -11.1121 & 1.2015 & -0.3519 \tabularnewline
			0.739 & -47.83 & 1.07 & -12.3369 & 1.2763 & -0.4562 \tabularnewline
			0.739 & -47.82 & 1.06 & -13.3770 & 1.3399 & -0.5406 \tabularnewline \hline
		\end{tabular}
	\end{center}
	\caption{$0.739$ [Nm] nagyságú nyomatékhoz tartozó adatok}
	\label{tablazat2}
\end{table}

A táblázatok a nyomaték nagysága, a kikapcsolási szög és a fellendítési idő mellett tartalmazzák a belső inga felső, függőleges helyzete ($\mathit{\Phi}_1=0\ [\degree]$) esetén fennálló inga-elemek közti szögkülönbséget, valamint a szögsebességeket. Megvizsgálva a kapott értékeket, láthatjuk, hogy több ideális paraméterkombinációt is találtunk. Azonban a szögkülönbséget és a szögsebességeket tekintve legkedvezőbb a \ref{tablazat1} táblázaton belüli, vastag betűvel jelzett eset. Tehát a belső csuklópontra egy $0.738$ [Nm] nagyságú nyomatékot kapcsolva, majd a belső inga $132\ [\degree]$-os, óramutató járásával ellentétes irányú elfordulása ($\mathit{\Phi}_1=-48\ [\degree]$) után lekapcsolva azt, elérhetjük, hogy a kettős inga a támasztott követelményeknek eleget téve, a kívánt módon lendüljön fel. Vagyis sikerült a lengőrendszert úgy a felső, instabil egyensúlyi helyzet kis környezetébe juttatni, hogy azt a későbbiekben egy lineáris szabályozó stabilizálni tudja.\par
A \ref{eq1} és \ref{eq2} nemlineáris differenciálegyenletekkel leírt, nyugalomban lévő ($\mathit{\Phi}_1=\mathit{\Phi}_2=-180\ [\degree]$, $\mathit{\omega_1}=\mathit{\omega_2}=0$ [rad/s]) mechanikai rendszert a kiválasztott nyomatékkal, a meghatározott pontig gerjesztve az alábbi eredményeket kapjuk.\\ [0.25 cm]
A kettős inga szögelfordulás függvényei:
\begin{figure}[h]
	\centering
	\subfloat[A belső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.465\linewidth]{fellenditesfi1}}} 
	\qquad
	\subfloat[A külső inga-elem szögelfordulás függvénye]{{\includegraphics[width=0.465\linewidth]{fellenditesfi2}}}
	\caption{A fellendítés során kapott szögelfordulás függvények}
	\label{fellenditesszog}
\end{figure}\\
A kettős inga szögsebesség függvényei:
\begin{figure}[h]
	\centering
	\subfloat[A belső inga-elem szögsebesség függvénye]{{\includegraphics[width=0.465\linewidth]{fellenditesomega1}}} 
	\qquad
	\subfloat[A külső inga-elem szögsebesség függvénye]{{\includegraphics[width=0.465\linewidth]{fellenditesomega2}}}
	\caption{A fellendítés során kapott szögsebesség függvények}
	\label{fellenditesszogseb}
\end{figure}\\
A trajektória vetületei: \vspace{-5 mm}
\begin{figure}[h]
	\centering
	\subfloat[A trajektória vetülete a szögelfordulások síkjára]{{\includegraphics[width=0.505\linewidth]{szogelftraj}}} 
	\qquad\qquad
	\subfloat[A trajektória vetülete a szögsebességek síkjára]{{\includegraphics[width=0.28\linewidth]{szogsebtraj}}}
	\caption{A fellendítés trajektóriájának vetületei}
	\label{trajektoriak}
\end{figure}\\
A szögelfordulás-, valamint a szögsebesség függvények ábráján látható első függőleges vonal a nyomaték kikapcsolásának időpontját jelzi, míg a második azt a pontot mutatja, amikor a belső inga eléri a felső, függőleges helyzetet. Megfigyelhető, hogy a nyomaték kikapcsolása után a szögsebességek fokozatosan csökkennek, s a két inga-elem minimális szögkülönbséggel, együtt halad a felső pozíció felé. Mind a \ref{fellenditesszog} és \ref{fellenditesszogseb} ábra, mind a trajektória vetületei szemléltetik, hogy a rendszer az instabil egyensúlyi helyzet környezetében meglehetősen kis szögsebességekkel rendelkezik, továbbá az elemek közti szögeltérés is optimális.\\ [0.25 cm]
Az alkalmazott gerjesztőnyomaték függvénye:
\begin{figure}[h]
	\centering
	\includegraphics[width=0.545\linewidth]{fellenditesinyomatek}
	\caption{A fellendítés nyomatékának függvénye}
	\label{fellenditonyomatek}
\end{figure}\par
A nyomatéknak és a kikapcsolási szögnek a szimuláció, hangolás során talált értéke elsősorban a rendszer mechanikai, illetve geometriai paramétereitől függ. A módszer azonban más kettős ingák fellendítésére ugyanúgy alkalmas lehet.\par
A fejezet zárásaként megtekinthetjük a lengőrendszer fellendítésének folyamatát a lenti képsorozaton is. A kettős inga $t=0$ [s]-ban nyugalomból indul. A mozgást $t=1.23$ [s]-ig, azaz $\mathit{\Phi}_1=0\  [\degree]$-ig követhetjük azonos, $\Delta t=0.15375$ [s] nagyságú időlépésekkel. \vspace{2 mm}

\begin{figure}[!htb]
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{1}
	\caption*{1. $t=0$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{2}
	\caption*{2. $t=0.15375$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}%
	\includegraphics[width=0.97\linewidth]{3}
	\caption*{3. $t=0.3075$ [s]}
	\endminipage \\ [0.3 cm]
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{4}
	\caption*{4. $t=0.46125$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{5}
	\caption*{5. $t=0.615$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}%
	\includegraphics[width=0.97\linewidth]{6}
	\caption*{6. $t=0.76875$ [s]}
	\endminipage \\ [0.3 cm]
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{7}
	\caption*{7. $t=0.9225$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}
	\includegraphics[width=0.97\linewidth]{8}
	\caption*{8. $t=1.07625$ [s]}
	\endminipage\hfill
	\minipage{0.32\textwidth}%
	\includegraphics[width=0.97\linewidth]{9}
	\caption*{9. $t=1.23$ [s]}
	\endminipage \\ [0.3 cm]
	\caption{A kettős inga fellendítésének folyamata a t $\in\, [0,1.23]$ [s] időtartományon, $\Delta t=0.15375$ [s]-os időlépésekkel}
\end{figure}

\end{section}

\chapter{A kettős inga szabályozása}
%PÁRATLAN OLDALON!

\lhead{ \textit{5. Fejezet - A kettős inga szabályozása}}
\rhead{ \textit{ \thepage}}

Ebben a fejezetben

%Először pólusok -> instabil jelleg -> szabályozzuk. szabályozható? (irányíthatóság, megfigyelhetőség)

%Mintavételi idő! Shannon-féle mintavételi törvény!

\end{document}
