!/
!/ Author : Toru Shiozaki
!/ Machine generated code
!/
      subroutine breitroot1(ta, rr, ww, n)
      implicit none
      integer i, j, n, offset, it, boxof
      double precision t, t2, d, e, f, g
      double precision ta(*), rr(*), ww(*)
      double precision ax(1)
      double precision aw(1)
      data (ax(i), i = 1, 1)/
     &   1.500000000000000d+00
     &/
      data (aw(i), i = 1, 1)/
     &   4.431134627263790d-01
     &/
      double precision x(384)
      double precision w(384)
      data x /
     &  1.057986693394197d+00,
     & -7.144282222753716d-02,
     & -2.054859453327102d-04,
     &  2.331146756779162d-04,
     &  1.331077830017425d-06,
     & -1.115171509111191d-06,
     & -8.999527162393619d-09,
     &  5.640469957441268d-09,
     &  5.868492293682324d-11,
     & -2.886120987688697d-11,
     & -3.681089987315752d-13,
     &  1.487463219474723d-13,
     &  7.851220690400758d-01,
     & -6.306370112940925d-02,
     &  2.077090515664856d-03,
     &  1.107583656237493d-04,
     & -1.277451547927649d-05,
     & -1.044487389335647d-08,
     &  6.009838083002671d-08,
     & -2.175115414441252d-09,
     & -2.080059059032380d-10,
     &  1.804354748413915d-11,
     &  3.447550363921192d-13,
     & -9.980335735361091d-14,
     &  5.700728455141461d-01,
     & -4.424389871474052d-02,
     &  2.353767967817790d-03,
     & -4.387325740927295d-05,
     & -5.008798807766955d-06,
     &  4.884288295089779d-07,
     & -1.173290746093463d-08,
     & -1.185932002907759d-09,
     &  1.202903061066997d-10,
     & -2.735219425546217d-12,
     & -3.083533159109398d-13,
     &  2.983920551688735d-14,
     &  4.263165133101116d-01,
     & -2.826877076462550d-02,
     &  1.603798639833974d-03,
     & -6.621801200235662d-05,
     &  9.705629527727786d-07,
     &  1.168434546975322d-07,
     & -1.170416650328710d-08,
     &  5.063515401445807d-10,
     & -9.010812787342137d-13,
     & -1.489830494944563d-12,
     &  1.095625215815169d-13,
     & -3.294304605340185d-15,
     &  3.343787695313898d-01,
     & -1.825217046123134d-02,
     &  9.421690541285707d-04,
     & -4.289019050586966d-05,
     &  1.506333938826640d-06,
     & -2.375199158095306d-08,
     & -1.642802891862454d-09,
     &  1.802845072073252d-10,
     & -9.557839910687163d-12,
     &  2.790196517090680d-13,
     &  1.998323378968711d-15,
     & -7.883517764092141d-16,
     &  2.737062959617326d-01,
     & -1.240744574659713d-02,
     &  5.529249190801726d-04,
     & -2.352162653738726d-05,
     &  9.051198513379506d-07,
     & -2.839050037196454d-08,
     &  5.109421933983436d-10,
     &  1.407039538947661d-11,
     & -1.949673343458228d-12,
     &  1.127750829710148d-13,
     & -4.374187315635119d-15,
     &  1.015963850541615d-16,
     &  2.314321486648044d-01,
     & -8.905431552066596d-03,
     &  3.411221934963245d-04,
     & -1.287092914971452d-05,
     &  4.676183879261686d-07,
     & -1.569514089798894d-08,
     &  4.499257284264834d-10,
     & -8.827025894078617d-12,
     & -4.528147109906610d-14,
     &  1.575661232617235d-14,
     & -9.771178325583664d-16,
     &  4.140295847801633d-17,
     &  2.004425977723634d-01,
     & -6.687523573427587d-03,
     &  2.228778169247477d-04,
     & -7.395888424548810d-06,
     &  2.423064077452070d-07,
     & -7.701382214546719d-09,
     &  2.300343695813390d-10,
     & -6.089986328575896d-12,
     &  1.239443061683193d-13,
     & -7.602553606340901d-16,
     & -9.582281811754662d-17,
     &  6.715820944592601d-18,
     &  1.767762141109114d-01,
     & -5.203610350329985d-03,
     &  1.531372602004397d-04,
     & -4.501650850226160d-06,
     &  1.318227678779673d-07,
     & -3.819713632881610d-09,
     &  1.080449216432170d-10,
     & -2.911511025958435d-12,
     &  7.155059166417748d-14,
     & -1.459122203961374d-15,
     &  1.714723604801397d-17,
     &  3.878715586059126d-19,
     &  1.581138140629937d-01,
     & -4.163747036495255d-03,
     &  1.096420218213107d-04,
     & -2.886384528312645d-06,
     &  7.590574699775137d-08,
     & -1.989587919128444d-09,
     &  5.170552283026493d-11,
     & -1.318402418151381d-12,
     &  3.237169261565035d-14,
     & -7.407527938951930d-16,
     &  1.480711182046233d-17,
     & -2.146133793937241d-19,
     &  1.430193785658226d-01,
     & -3.407151823652130d-03,
     &  8.116779623216395d-05,
     & -1.933526852665239d-06,
     &  4.604702000492541d-08,
     & -1.095581820289707d-09,
     &  2.599519200175242d-11,
     & -6.125793140220651d-13,
     &  1.422136136841771d-14,
     & -3.206123921352768d-16,
     &  6.847741816300389d-18,
     & -1.323580050468122d-19,
     &  1.305582405744536d-01,
     & -2.839564625134561d-03,
     &  6.175873588666825d-05,
     & -1.343196519956186d-06,
     &  2.921147589499571d-08,
     & -6.351259069163048d-10,
     &  1.379793425287233d-11,
     & -2.990809553309502d-13,
     &  6.447556601376684d-15,
     & -1.373798642932894d-16,
     &  2.861296678012272d-18,
     & -5.715160818258370d-20,
     &  1.200961151572054d-01,
     & -2.402883759028992d-03,
     &  4.807689503210470d-05,
     & -9.619199787150217d-07,
     &  1.924577355659989d-08,
     & -3.850393897402355d-10,
     &  7.701558903062076d-12,
     & -1.539414978130730d-13,
     &  3.071428657261740d-15,
     & -6.101756808188811d-17,
     &  1.201180640067729d-18,
     & -2.322474804736961d-20,
     &  1.111873974714942d-01,
     & -2.059732465291284d-03,
     &  3.815628079219078d-05,
     & -7.068398697088741d-07,
     &  1.309407210509106d-08,
     & -2.425616628010063d-10,
     &  4.493087025526283d-12,
     & -8.321161507172306d-14,
     &  1.540203208176981d-15,
     & -2.846673371464526d-17,
     &  5.243582978061400d-19,
     & -9.587378206994651d-21,
     &  1.035098338974698d-01,
     & -1.785183137531318d-03,
     &  3.078817421220887d-05,
     & -5.309884233881554d-07,
     &  9.157689161265787d-09,
     & -1.579375283561523d-10,
     &  2.723822044721043d-12,
     & -4.697317805182187d-14,
     &  8.099353198328178d-16,
     & -1.395889191299331d-17,
     &  2.402965421906549d-19,
     & -4.124487881562381d-21,
     &  9.682458365464176d-02,
     & -1.562093310508541d-03,
     &  2.520161113926371d-05,
     & -4.065833892823905d-07,
     &  6.559502515141727d-09,
     & -1.058258753430705d-10,
     &  1.707306223526556d-12,
     & -2.754389672061436d-14,
     &  4.443447182854592d-16,
     & -7.167301793964512d-18,
     &  1.155660357036138d-19,
     & -1.861209421046733d-21,
     &  9.095085938854892d-02,
     & -1.378359824587057d-03,
     &  2.088903632332902d-05,
     & -3.165732411070342d-07,
     &  4.797665769217359d-09,
     & -7.270858827993690d-11,
     &  1.101897296636707d-12,
     & -1.669918096052681d-14,
     &  2.530720658696433d-16,
     & -3.835100826358939d-18,
     &  5.811134365390768d-20,
     & -8.800694329041991d-22,
     &  8.574929257124383d-02,
     & -1.225239993899981d-03,
     &  1.750700207026308d-05,
     & -2.501510910788369d-07,
     &  3.574316597596288d-09,
     & -5.107208887892981d-11,
     &  7.297500866742994d-13,
     & -1.042712007914412d-14,
     &  1.489887058346579d-16,
     & -2.128815092662221d-18,
     &  3.041646382120743d-20,
     & -4.344618678335355d-22,
     &  8.111071056537980d-02,
     & -1.096290919106305d-03,
     &  1.481744853332377d-05,
     & -2.002723703977008d-07,
     &  2.706877787982895d-09,
     & -3.658611170781422d-11,
     &  4.944972139946927d-13,
     & -6.683614174622216d-15,
     &  9.033552960658669d-17,
     & -1.220969179453738d-18,
     &  1.650239748312328d-20,
     & -2.229968378756115d-22,
     &  7.694837640638635d-02,
     & -9.866798490756386d-04,
     &  1.265182152029981d-05,
     & -1.622295093272588d-07,
     &  2.080207474479467d-09,
     & -2.667371152017094d-11,
     &  3.420268843585154d-13,
     & -4.385680836693997d-15,
     &  5.623591072936953d-17,
     & -7.210911539448567d-19,
     &  9.246249596817101d-21,
     & -1.185404358058847d-22,
     &  7.319250547113997d-02,
     & -8.927243167395140d-04,
     &  1.088850150118277d-05,
     & -1.328063577050375d-07,
     &  1.619830666731022d-09,
     & -1.975697122970543d-11,
     &  2.409745164363792d-13,
     & -2.939150766386149d-15,
     &  3.584863264978267d-17,
     & -4.372433983921355d-19,
     &  5.333025633602838d-21,
     & -6.503672283796038d-23,
     &  6.978631577988531d-02,
     & -8.115785350685460d-04,
     &  9.438236009785576d-06,
     & -1.097617730474387d-07,
     &  1.276472299480924d-09,
     & -1.484470855413484d-11,
     &  1.726362351034134d-13,
     & -2.007669568920946d-15,
     &  2.334815204849949d-17,
     & -2.715268456353038d-19,
     &  3.157715275344596d-21,
     & -3.671759456771296d-23,
     &  6.668313367115805d-02,
     & -7.410152021126283d-04,
     &  8.234518978515185d-06,
     & -9.150595374319278d-08,
     &  1.016858372940869d-09,
     & -1.129982157788204d-11,
     &  1.255690773493143d-13,
     & -1.395384260947041d-15,
     &  1.550618411413408d-17,
     & -1.723122086361966d-19,
     &  1.914816411580889d-21,
     & -2.127573437600437d-23,
     &  6.384423980690615d-02,
     & -6.792709245887912d-04,
     &  7.227104440231730d-06,
     & -7.689279299218559d-08,
     &  8.181010338282292d-10,
     & -8.704187681385158d-12,
     &  9.260822326244377d-14,
     & -9.853053873845546d-16,
     &  1.048315875054867d-17,
     & -1.115355894803848d-19,
     &  1.186683113722182d-21,
     & -1.262428794520145d-23,
     &  6.123724356957945d-02,
     & -6.249349093931703d-04,
     &  6.377550951236128d-06,
     & -6.508382797035222d-08,
     &  6.641898584053275d-10,
     & -6.778153371822801d-12,
     &  6.917203349381316d-14,
     & -7.059105858388098d-16,
     &  7.203919416475140d-18,
     & -7.351703739951724d-20,
     &  7.502519680997237d-22,
     & -7.655632257320562d-24,
     &  5.883484054145521d-02,
     & -5.768676142156884d-04,
     &  5.656108545011271d-06,
     & -5.545737546117121d-08,
     &  5.437520281950621d-10,
     & -5.331414725409118d-12,
     &  5.227379669488312d-14,
     & -5.125374711271604d-16,
     &  5.025360236198257d-18,
     & -4.927297402477483d-20,
     &  4.831148080528969d-22,
     & -4.736419676927769d-24,
     &  5.661385170722979d-02,
     & -5.341404831787391d-04,
     &  5.039509716559036d-06,
     & -4.754677651870560d-08,
     &  4.485944237574224d-10,
     & -4.232399581222562d-12,
     &  3.993185217304065d-14,
     & -3.767491200602163d-16,
     &  3.554553363832264d-18,
     & -3.353650730252294d-20,
     &  3.164103047362650d-22,
     & -2.985002873487706d-24,
     &  5.455447255899810d-02,
     & -4.959907448952424d-04,
     &  4.509379478569656d-06,
     & -4.099774741974261d-08,
     &  3.727376020317001d-10,
     & -3.388803744408602d-12,
     &  3.080985324668412d-14,
     & -2.801127266954748d-16,
     &  2.546689821223788d-18,
     & -2.315363932946449d-20,
     &  2.105050273731409d-22,
     & -1.913682066835144d-24,
     &  5.263968047576973d-02,
     & -4.617871188746624d-04,
     &  4.051075941783504d-06,
     & -3.553848865704591d-08,
     &  3.117651197303763d-10,
     & -2.734992216986845d-12,
     &  2.399300612412218d-14,
     & -2.104811630895123d-16,
     &  1.846468082670482d-18,
     & -1.619833494965560d-20,
     &  1.421015917085596d-22,
     & -1.246505190398957d-24,
     &  5.085476277156078d-02,
     & -4.310035220859456d-04,
     &  3.652834580803000d-06,
     & -3.095844880833592d-08,
     &  2.623785806385094d-10,
     & -2.223707008192934d-12,
     &  1.884632825687528d-14,
     & -1.597261183497956d-16,
     &  1.353708400668723d-18,
     & -1.147292911749707d-20,
     &  9.723519614044817d-23,
     & -8.240270859322372d-25,
     &  4.918693768379647d-02,
     & -4.031987115847757d-04,
     &  3.305129546155459d-06,
     & -2.709304619038437d-08,
     &  2.220890714338660d-10,
     & -1.820524547286319d-12,
     &  1.492333506495387d-14,
     & -1.223306380531084d-16,
     &  1.002777525355174d-18,
     & -8.220040223452146d-21,
     &  6.738190601767012d-23,
     & -5.523107050623459d-25,
     &  4.762504762507144d-02,
     & -3.780003795005700d-04,
     &  3.000192000382501d-06,
     & -2.381254762508930d-08,
     &  1.890003787507110d-10,
     & -1.500097500289126d-12,
     &  1.190628571883334d-14,
     & -9.450028387566298d-17,
     &  7.500495001942498d-19,
     & -5.953148812565445d-21,
     &  4.725018900051778d-23,
     & -3.750015000040365d-25
     &/
      data w /
     &  4.135653584758195d-01,
     & -1.066439704592393d-01,
     &  1.752146780833233d-02,
     & -2.154035108083543d-03,
     &  2.123679824895491d-04,
     & -1.749369492647350d-05,
     &  1.237819265129784d-06,
     & -7.676522952021522d-08,
     &  4.237208259260642d-09,
     & -2.107052960597916d-10,
     &  9.532877171342738d-12,
     & -3.950665544868368d-13,
     &  1.599461505363333d-01,
     & -3.096198038212439d-02,
     &  4.255772773271153d-03,
     & -4.654909238812563d-04,
     &  4.233489892280168d-05,
     & -3.289699797763587d-06,
     &  2.227993127460179d-07,
     & -1.335652563041112d-08,
     &  7.176391768087074d-10,
     & -3.491410473430721d-11,
     &  1.551298009286964d-12,
     & -6.332303946483246d-14,
     &  8.026264548904928d-02,
     & -1.133995737427459d-02,
     &  1.243259373351908d-03,
     & -1.159576319402131d-04,
     &  9.415673853935415d-06,
     & -6.737652179339188d-07,
     &  4.291072677258729d-08,
     & -2.454290753729082d-09,
     &  1.271030018274720d-10,
     & -6.004510636076742d-12,
     &  2.604804363995930d-13,
     & -1.042458350064328d-14,
     &  4.859569796092199d-02,
     & -5.148911980098271d-03,
     &  4.438751383843736d-04,
     & -3.424488814882181d-05,
     &  2.404895705877444d-06,
     & -1.541735418447310d-07,
     &  9.032246467106875d-09,
     & -4.844980385243089d-10,
     &  2.386691820538419d-11,
     & -1.083687867755238d-12,
     &  4.553481216480986d-14,
     & -1.775509251352928d-15,
     &  3.318841541503713d-02,
     & -2.763244316554662d-03,
     &  1.905338480625239d-04,
     & -1.207874145926498d-05,
     &  7.193505769630907d-07,
     & -4.034421045107958d-08,
     &  2.125092496484539d-09,
     & -1.048140678382111d-10,
     &  4.831960683319789d-12,
     & -2.081263277623661d-13,
     &  8.382481202901493d-15,
     & -3.158203880028124d-16,
     &  2.447924039107278d-02,
     & -1.670337154526220d-03,
     &  9.487009895039248d-05,
     & -5.006890236939414d-06,
     &  2.523646848363528d-07,
     & -1.222934372643786d-08,
     &  5.691127493146268d-10,
     & -2.533371448326467d-11,
     &  1.074109554554316d-12,
     & -4.322752799410106d-14,
     &  1.647705940116275d-15,
     & -5.935895690276484d-17,
     &  1.901257983093334d-02,
     & -1.097797176431922d-03,
     &  5.282346462942369d-05,
     & -2.369704966560619d-06,
     &  1.022035381291580d-07,
     & -4.282846433782955d-09,
     &  1.747862721399167d-10,
     & -6.934187011465713d-12,
     &  2.664053821042440d-13,
     & -9.868867536559238d-15,
     &  3.511211871133835d-16,
     & -1.194910648313331d-17,
     &  1.531866754435402d-02,
     & -7.664545839261231d-04,
     &  3.196393643077859d-05,
     & -1.243885426763538d-06,
     &  4.663965470028716d-08,
     & -1.706421432978073d-09,
     &  6.122804189144124d-11,
     & -2.156405613094402d-12,
     &  7.441977150769849d-14,
     & -2.508524188712748d-15,
     &  8.226685739996370d-17,
     & -2.612347731342939d-18,
     &  1.268477574194269d-02,
     & -5.599240346078043d-04,
     &  2.060071644229961d-05,
     & -7.074073551965400d-07,
     &  2.341969436362825d-08,
     & -7.576730424366755d-10,
     &  2.410717114225890d-11,
     & -7.564247508605618d-13,
     &  2.341539700006752d-14,
     & -7.141043632992439d-16,
     &  2.140034617675300d-17,
     & -6.276204210715181d-19,
     &  1.072862287662884d-02,
     & -4.236815957127729d-04,
     &  1.394531016987716d-05,
     & -4.284172633161849d-07,
     &  1.269106180327052d-08,
     & -3.675409696214260d-10,
     &  1.047863069206222d-11,
     & -2.951796290593520d-13,
     &  8.229608960255464d-15,
     & -2.271270305495372d-16,
     &  6.198869790374740d-18,
     & -1.668561465030425d-19,
     &  9.228695623565627d-03,
     & -3.297131299395551d-04,
     &  9.817758366971190d-06,
     & -2.728598017704787d-07,
     &  7.312637893641860d-09,
     & -1.916176195022723d-10,
     &  4.944459162433632d-12,
     & -1.261469327304046d-13,
     &  3.189406294216692d-15,
     & -8.000255819943310d-17,
     &  1.991234616332444d-18,
     & -4.911271281247086d-20,
     &  8.048660322776655d-03,
     & -2.625338887107062d-04,
     &  7.137025038788259d-06,
     & -1.810916752028764d-07,
     &  4.430876648208031d-09,
     & -1.060033991982483d-10,
     &  2.497516150226103d-12,
     & -5.819196113628237d-14,
     &  1.344291859419495d-15,
     & -3.083708903932189d-17,
     &  7.029905834445299d-19,
     & -1.592076847859000d-20,
     &  7.100466902341420d-03,
     & -2.130672922355064d-04,
     &  5.328547921968880d-06,
     & -1.243794327640887d-07,
     &  2.799613986105193d-09,
     & -6.161539611874477d-11,
     &  1.335512612964982d-12,
     & -2.862854052218697d-14,
     &  6.085420723676045d-16,
     & -1.284902965904507d-17,
     &  2.697837302973300d-19,
     & -5.633808751437497d-21,
     &  6.324966011972775d-03,
     & -1.757311768923202d-04,
     &  4.069072701187250d-06,
     & -8.794034532608399d-08,
     &  1.832695466957234d-09,
     & -3.734515753007592d-11,
     &  7.494593519896410d-13,
     & -1.487514839958279d-14,
     &  2.927747683813691d-16,
     & -5.724526544293141d-18,
     &  1.113289752186427d-19,
     & -2.154426347039960d-21,
     &  5.681099140852745d-03,
     & -1.469522872531116d-04,
     &  3.167899377741413d-06,
     & -6.373988222908077d-08,
     &  1.236687839418025d-09,
     & -2.346122399237596d-11,
     &  4.383403336434155d-13,
     & -8.099782150469673d-15,
     &  1.484224295599604d-16,
     & -2.701922771248382d-18,
     &  4.892614674764565d-20,
     & -8.817653407041110d-22,
     &  5.139564276280682d-03,
     & -1.243645221694533d-04,
     &  2.507920247977568d-06,
     & -4.720356285170548d-08,
     &  8.567308298362593d-10,
     & -1.520392014730100d-11,
     &  2.657277809894234d-13,
     & -4.593247384237508d-15,
     &  7.873508869302429d-17,
     & -1.340811796480891d-18,
     &  2.271285871687245d-20,
     & -3.829655292905190d-22,
     &  4.678946550556683d-03,
     & -1.063549572339188d-04,
     &  2.014703217768674d-06,
     & -3.562110879848837d-08,
     &  6.073122120068501d-10,
     & -1.012413675034798d-11,
     &  1.662164942930265d-13,
     & -2.698930916722854d-15,
     &  4.345855943362673d-17,
     & -6.952018235722285d-19,
     &  1.106251230120596d-20,
     & -1.752273553681212d-22,
     &  4.283269318842086d-03,
     & -9.179605331848567d-05,
     &  1.639508063160909d-06,
     & -2.733036443695931d-08,
     &  4.393243128554762d-10,
     & -6.905042006491795d-12,
     &  1.068851730258372d-13,
     & -1.636327055875213d-15,
     &  2.484212332861477d-17,
     & -3.746791908434509d-19,
     &  5.621322243996755d-21,
     & -8.395309013619787d-23,
     &  3.940396043305408d-03,
     & -7.988201151284797d-05,
     &  1.349573985466081d-06,
     & -2.128071127094307d-08,
     &  3.235812117714724d-10,
     & -4.810847708703187d-12,
     &  7.044171304931821d-14,
     & -1.020092791670516d-15,
     &  1.464925013993331d-17,
     & -2.089986176400013d-19,
     &  2.966058226349572d-21,
     & -4.190307323916755d-23,
     &  3.640959158449784d-03,
     & -7.002564001642558d-05,
     &  1.122367213065414d-06,
     & -1.679013143131793d-08,
     &  2.422036803923001d-10,
     & -3.416240491374276d-12,
     &  4.745544659255305d-14,
     & -6.519662724509052d-16,
     &  8.882395784650675d-18,
     & -1.202229591640969d-19,
     &  1.618650135487074d-21,
     & -2.169487867635906d-23,
     &  3.377625055343649d-03,
     & -6.179156849562061d-05,
     &  9.420672508340389d-07,
     & -1.340526515447216d-08,
     &  1.839400324895536d-10,
     & -2.467846138017629d-12,
     &  3.260842198355359d-14,
     & -4.261307664627624d-16,
     &  5.522324302321272d-18,
     & -7.109731733863023d-20,
     &  9.105269551140079d-22,
     & -1.160856254088352d-23,
     &  3.144579069692897d-03,
     & -5.485194517482078d-05,
     &  7.973609904840891d-07,
     & -1.081828414100072d-08,
     &  1.415366735605166d-10,
     & -1.810591519744941d-12,
     &  2.281086845722602d-14,
     & -2.842264806526018d-16,
     &  3.511988985652421d-18,
     & -4.311158917850439d-20,
     &  5.264331156728034d-22,
     & -6.399481196589546d-24,
     &  2.937157824996331d-03,
     & -4.895640832540258d-05,
     &  6.800235834725452d-07,
     & -8.816141034468411d-09,
     &  1.102148610266473d-10,
     & -1.347232710227029d-12,
     &  1.621865681798140d-14,
     & -1.931027712709452d-16,
     &  2.279964073211056d-18,
     & -2.674359393942946d-20,
     &  3.120468565482904d-22,
     & -3.624747587604266d-24,
     &  2.751582327594653d-03,
     & -4.391133493588563d-05,
     &  5.839851625408772d-07,
     & -7.248811980568326d-09,
     &  8.676385087918351d-11,
     & -1.015433718747613d-12,
     &  1.170399463984184d-14,
     & -1.334190931564087d-16,
     &  1.508230490745103d-18,
     & -1.693829560234687d-20,
     &  1.892256345942885d-22,
     & -2.104519600504218d-24,
     &  2.584761240793730d-03,
     & -3.956524707720406d-05,
     &  5.047047508050445d-07,
     & -6.008976519042477d-09,
     &  6.898751233313386d-11,
     & -7.744282563821105d-13,
     &  8.561732222065785d-15,
     & -9.361456201905489d-17,
     &  1.015058398462874d-18,
     & -1.093429783498432d-20,
     &  1.171653073097142d-22,
     & -1.249896817374903d-24,
     &  2.434143871810106d-03,
     & -3.579838410794356d-05,
     &  4.387425901699890d-07,
     & -5.018749879238761d-09,
     &  5.535898083967781d-11,
     & -5.970645754649600d-13,
     &  6.341970867074628d-15,
     & -6.662365953662773d-17,
     &  6.940624760109292d-19,
     & -7.183248492750280d-21,
     &  7.395226454323649d-23,
     & -7.579713870401952d-25,
     &  2.297608917210032d-03,
     & -3.251514252254305d-05,
     &  3.834631490292616d-07,
     & -4.220858563891419d-09,
     &  4.480068602231937d-11,
     & -4.649531291505860d-13,
     &  4.752294900901216d-15,
     & -4.803954443047130d-17,
     &  4.815708693534111d-19,
     & -4.795940360253550d-21,
     &  4.751115597497762d-23,
     & -4.685883728220409d-25,
     &  2.173379292042435d-03,
     & -2.963852137211661d-05,
     &  3.368257383492216d-07,
     & -3.572671024824223d-09,
     &  3.654158788717742d-11,
     & -3.654453283032621d-13,
     &  3.599374959628388d-15,
     & -3.506170609690566d-17,
     &  3.386919083815422d-19,
     & -3.250340274810860d-21,
     &  3.102852216804809d-23,
     & -2.948970470639993d-25,
     &  2.059956231480111d-03,
     & -2.710599091783694d-05,
     &  2.972348264873238d-07,
     & -3.042096344982233d-09,
     &  3.002291131612693d-11,
     & -2.897164953813359d-13,
     &  2.753361925182086d-15,
     & -2.587942652059519d-17,
     &  2.412191645113400d-19,
     & -2.233681237769136d-21,
     &  2.057495218055001d-23,
     & -1.886844963996866d-25,
     &  1.956067801881088d-03,
     & -2.486638490378408d-05,
     &  2.634316497790524d-07,
     & -2.604725597432706d-09,
     &  2.483490579790712d-11,
     & -2.315280471372690d-13,
     &  2.125760555459292d-15,
     & -1.930306780552970d-17,
     &  1.738217573602923d-19,
     & -1.555013481807576d-21,
     &  1.383797140440405d-23,
     & -1.226006718421754d-25,
     &  1.860628303764434d-03,
     & -2.287753822164470d-05,
     &  2.344147880046653d-07,
     & -2.241812684402721d-09,
     &  2.067379002693027d-11,
     & -1.864152370643060d-13,
     &  1.655435802264321d-15,
     & -1.453933445336660d-17,
     &  1.266317219503761d-19,
     & -1.095702660121998d-21,
     &  9.430857108371231d-24,
     & -8.081536494792655d-26,
     &  1.772705991748090d-03,
     & -2.110447363715505d-05,
     &  2.093813188035964d-07,
     & -1.938830412437076d-09,
     &  1.731203545102832d-11,
     & -1.511461004080499d-13,
     &  1.299616981265166d-15,
     & -1.105185166202646d-17,
     &  9.320098747539451d-20,
     & -7.808330630379991d-22,
     &  6.507349278507902d-24,
     & -5.399286049308797d-26
     &/
      offset = -1
      do i=1, n
        t = ta(i)
        offset = offset + 1
        if (t < 0.0d0) then
          rr(offset+1:offset+1) = 0.5d0
          ww(offset+1:offset+1) = 0.0d0
        else if (t >= 64.0d0) then
          t = 1.0d0/dsqrt(t)
          rr(offset+1:offset+1) = ax(1:1)*t*t
          ww(offset+1:offset+1) = aw(1:1)*t*t*t
        else
          it = int(t*   0.500000000000000d0)
          t = (t-it*2.000000000000000-   1.000000000000000d0)
     &     *   1.000000000000000d0
          t2 = t * 2.0d0
          do j=1, 1
            boxof = it*12+12*(j-1)
            d = x(boxof+12)
            e = w(boxof+12)
            f = t2*d + x(boxof+11)
            g = t2*e + w(boxof+11)
            d = t2*f - d + x(boxof+10)
            e = t2*g - e + w(boxof+10)
            f = t2*d - f + x(boxof+9)
            g = t2*e - g + w(boxof+9)
            d = t2*f - d + x(boxof+8)
            e = t2*g - e + w(boxof+8)
            f = t2*d - f + x(boxof+7)
            g = t2*e - g + w(boxof+7)
            d = t2*f - d + x(boxof+6)
            e = t2*g - e + w(boxof+6)
            f = t2*d - f + x(boxof+5)
            g = t2*e - g + w(boxof+5)
            d = t2*f - d + x(boxof+4)
            e = t2*g - e + w(boxof+4)
            f = t2*d - f + x(boxof+3)
            g = t2*e - g + w(boxof+3)
            d = t2*f - d + x(boxof+2)
            e = t2*g - e + w(boxof+2)
            rr(offset+j) = t*d - f + x(boxof+1)*0.5d0
            ww(offset+j) = t*e - g + w(boxof+1)*0.5d0
          enddo
        endif
      enddo
      end subroutine
