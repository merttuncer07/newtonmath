# Slice 9 doğrulama raporu

Tarih: 2026-10-01. Ortam: yerel macOS, Apple Clang 21.0.0, C11; `/usr/bin/bc`.
Başlangıç: kabul edilmiş Slice 8, `15e486d`. Rasyonel integrasyon, kısmi kesirler,
konik alan değerleri ve belirli integral hesabı tamamlandı. Sonuçlar aşağıda mevcut program çıktılarıyla gösteriliyor.

## Verilen 20 cevap

İlk sütun çalıştırılan ifade, ikinci sütun sağlanan cevap anahtarı, üçüncü sütun gerçek çıktı.
Sabit farkı ve terim sırası serbest; her ilkelin türevi program içinde tam aritmetikte denetleniyor.
`pisqrt(3)` mevcut yazıcıda `pi*sqrt(3)` çarpımıdır. Her satır `[exact]` verdi.

| İfade | Cevap anahtarı | Program çıktısı |
|---|---|---|
| `integral(1/(x^2 - 1), x)` | `log\|x - 1\|/2 - log\|x + 1\|/2` | `(-1/2)*log\|x + 1\| + (1/2)*log\|x - 1\|  [exact]` |
| `integral(x/(x^2 + 1), x)` | `log\|x^2 + 1\|/2` | `(1/2)*log\|x^2 + 1\|  [exact]` |
| `integral(1/(x^2 + 1), x)` | `atan(x)` | `atan(x)  [exact]` |
| `integral(1/(x^2 + x + 1), x)` | `(2/sqrt(3))*atan((2x + 1)/sqrt(3))` | `(2sqrt(3)/3)*atan(2sqrt(3)x/3 + sqrt(3)/3)  [exact]` |
| `integral((x^3 + 1)/(x^2 - 1), x)` | `x^2/2 + log\|x - 1\|` | `x^2/2 + log\|x - 1\|  [exact]` |
| `integral(1/(x - 1)^2, x)` | `-1/(x - 1)` | `-1/(x - 1)  [exact]` |
| `integral(1/(x*(x + 1)^2), x)` | `log\|x\| - log\|x + 1\| + 1/(x + 1)` | `1/(x + 1) + (-1)*log\|x + 1\| + log\|x\|  [exact]` |
| `integral((2x + 3)/(x^2 + 2x + 5), x)` | `log\|x^2 + 2x + 5\| + atan((x + 1)/2)/2` | `log\|x^2 + 2x + 5\| + (1/2)*atan(x/2 + 1/2)  [exact]` |
| `integral(1/(x^2 - 2), x)` | `(sqrt(2)/4)*(log\|x - sqrt(2)\| - log\|x + sqrt(2)\|)` | `(sqrt(2)/4)*log\|-sqrt(2) + x\| + (-sqrt(2)/4)*log\|sqrt(2) + x\|  [exact]` |
| `integral(1/(x^4 - 1), x)` | `log\|x - 1\|/4 - log\|x + 1\|/4 - atan(x)/2` | `(-1/4)*log\|x + 1\| + (1/4)*log\|x - 1\| + (-1/2)*atan(x)  [exact]` |
| `integral(x^2/(x + a), x)` | `x^2/2 - a*x + a^2*log\|x + a\|` | `-ax + x^2/2 + (a^2)*log\|a + x\|  [exact]` |
| `integral(1/(x^3 + 1), x)` | `log\|x + 1\|/3 - log\|x^2 - x + 1\|/6 + atan((2x - 1)/sqrt(3))/sqrt(3)` | `(1/3)*log\|x + 1\| + (-1/6)*log\|x^2 - x + 1\| + (sqrt(3)/3)*atan(2sqrt(3)x/3 - sqrt(3)/3)  [exact]` |
| `apart(1/(x^4 - 1), x)` | `1/(4*(x - 1)) - 1/(4*(x + 1)) - 1/(2*(x^2 + 1))` | `(-1/4)/(x + 1) + (1/4)/(x - 1) + (-1/2)/(x^2 + 1)  [exact]` |
| `apart(1/(x*(x + 1)^2), x)` | `1/x - 1/(x + 1) - 1/(x + 1)^2` | `-1/(x + 1) + -1/(x^2 + 2x + 1) + 1/x  [exact]` |
| `apart(1/(x^3 + 1), x)` | `1/(3*(x + 1)) - (x - 2)/(3*(x^2 - x + 1))` | `(1/3)/(x + 1) + (-x/3 + 2/3)/(x^2 - x + 1)  [exact]` |
| `integral(1/(1 + x), x, 0, 1)` | `log(2)` | `log(2)  [exact]` |
| `integral(1/(1 + x^2), x, 0, 1)` | `pi/4` | `pi/4  [exact]` |
| `integral(1/(x^2 + x + 1), x, 0, 1)` | `sqrt(3)*pi/9` | `pisqrt(3)/9  [exact]` |
| `integral(x/(x^2 + 1), x, 0, 2)` | `log(5)/2` | `log(5)/2  [exact]` |
| `integral(1/(x^3 + 1), x, 0, 1)` | `log(2)/3 + sqrt(3)*pi/9` | `pisqrt(3)/9 + log(2)/3  [exact]` |

## İstenen retler

| İfade | Gerçek çıktı |
|---|---|
| `integral(1/(x^3 + x + 1), x)` | `error: remaining degree-3 factor has no supported linear/quadratic split; later: Rothstein-Trager` |
| `integral(1/(x^2 + a), x)` | `error: letters in a nonlinear conic factor: sign unknown; later` |
| `integral(1/x, x, -1, 1)` | `error: pole in the closed integration interval` |

Üçüncü derece kalan için program indirgenemezlik ispatı iddia etmiyor; desteklenen bir ayrışma bulamadığını
söylüyor. İkinci ret parametrik konik kısmın kapsam dışı olduğunu, üçüncüsü kapalı aralıktaki kutbu bildiriyor.

## 60 basamak

| Belirli integral | Programın yuvarlanmış sonucu |
|---|---|
| `1/(1+x), 0..1` | `0.693147180559945309417232121458176568075500134360255254120680` |
| `1/(1+x^2), 0..1` | `0.785398163397448309615660845819875721049292349843776455243736` |
| `1/(x^2+x+1), 0..1` | `0.604599788078072616864692752547385244094688749364246858523295` |
| `x/(x^2+1), 0..2` | `0.804718956217050187300379666613093819762800677134258860956324` |
| `1/(x^3+1), 0..1` | `0.835648848264721053337103459700110766786522127484331943230188` |

Beşi de `[bounded: 60 places guaranteed]` verdi ve bağımsız `bc -l` hesabıyla aynı 60 basamağa yuvarlandı.
Karşılaştırmada bc 75 basamakta çalışıyor; negatif sayılar da işaretine uygun yuvarlanıyor.

## Kontrol sayıları ve süre

- İzole integral modülü: **86 kontrol, 0 hata**. Bunların **30** tanesi türevle geri dönüş,
  **29** tanesi kısmi kesirleri toplayarak geri dönüş; ayrıca bağımsız beklenen ilkel biçimleri, parametre
  katsayıları, tekrarlar, büyük katsayılı çarpanlar, kutuplar ve kalıcılık denetleniyor.
- Dil fixture'ı: **109 incelenmiş çıktı**. İçinde **13** açık integral türevi farkı ve **5** açık apart toplamı
  farkı sıfırla karşılaştırılıyor. Böylece bu iki katmanda açıkça sayılan geri kontroller **43 türev**,
  **34 toplam**. Program ayrıca her kabul edilen integrali ve apart çağrısını kendi içinde kontrol ediyor;
  bu otomatik tekrarlar sayılara ikinci kez eklenmedi.
- **20 belirli integral** bağımsız bc ile 60 basamakta karşılaştırıldı: verilen beş cevap ve 15 ek örnek.
  Ek örnekler büyük/küçük log argümanlarını, genel atan değerlerini, ters/negatif/köklü sınırları,
  tekrarlı ikinci derece paydaları ve birden fazla ikinci derece çarpanı kapsıyor.
- **1 ek bc regresyonu**: kullanıcının `atan` adlı başka bir kuralı ve taşınmış serisi, matematiksel
  `atan(2)` sabitinin değerini değiştiremiyor. Toplam **21 yeni bc karşılaştırması**.
- Son `make test`: `run.sh` **126**, `check_z` **9065**; birim testler arith **19**, cplx **13**, elim **9**,
  integ **86**, linalg **13**, ratfun **788**, ratmat **553**. Hepsi geçti.
- Son derlenmiş durumdaki `make test`: **7,24 saniye** gerçek süre. Önceki tüm test ikililerini yeniden
  derleyen geçiş **19,26 saniye** sürdü (o geçişte integral birim testinde 84 kontrol vardı; son iki büyük
  katsayı kontrolü eklenip ayrıca derlendi ve son 7,24 saniyelik geçişte çalıştı). Her ikisi de 60 saniyenin altında.
- AddressSanitizer/UndefinedBehaviorSanitizer: 109 satırlık dil fixture'ının çıktısı birebir aynı, tanı yok.
  İzole modülün 84 kontrole kadarki sürümü de bu araçlarla temiz geçti. Son iki çarpan kontrolü normal
  derlemede geçti. Arena/kalıcı ayırma düzeni nedeniyle bu araçlarda leak raporlaması kapalıydı.
- `git diff --check`: temiz. Bu turda GNU GCC/Linux üzerinde yeniden çalıştırılmadı; ortam yukarıda yazılı.

## Eski .out satırları

**Değişen eski satır sayısı: 0.** Slice 8 başlangıcıyla karşılaştırmada `tests/series.out`, `newton.out`,
`irrational.out`, `rules.out`, `integration.out` ve `rational.out` dosyaları değişmedi. Bu nedenle önce/sonra
satırı ve değişiklik gerekçesi yok. Yeni `tests/integral.out` eklendi. Açık seri istekleri ve eski nonrational
integrasyon davranışı mevcut çıktıları koruyor.

## Yapılmayanlar ve sınırlar

- Rothstein-Trager ve genel yüksek dereceli cebirsel logaritma integrasyonu yok. Çarpan araması sınırlı;
  64'ten büyük payda derecesi, büyük rasyonel kök araması veya çok fazla bölen/adaya karşı açık ret var.
  Arama başarısızlığı indirgenemezlik ispatı değildir.
- Parametreli ikinci derece log/atan kısmı ve genel çok parametreli çarpanlama yok. Tek doğrusal çarpan,
  tekrarları ve x'ten bağımsız rasyonel parametre katsayıları destekleniyor.
- Konik alanların genel aritmetiği, yeniden integrali ve farklı bir parametreye göre türevi yok.
  Sayısal yerine koyma, kalan log/atan argümanlarının gerçek sayı olmasını gerektiriyor.
- Belirli integral sınırları kesin gerçek değerlerdir: rasyoneller, desteklenen kökler ve adlandırılmış konik
  sabitler. Doğrudan yaklaşık ball sınırları kabul edilmiyor. Köklü bir uç kutba çalışma hassasiyetinde
  ayrılamayacak kadar yakınsa sonuç tahmin edilmiyor.
- Genel log/atan özdeşlik sadeleştirmesi yok; tanıdık özel açılar pi'ye dönüşürken diğerleri `atan(k)` kalabilir.
- Rasyonel sadeleştirmede iptal edilen delikler korunmuyor; kapalı aralık denetimi sadeleşmiş fonksiyona ait.
- Kaynak notu: sağlanan tasarımın Problem IX/tablo atfı yereldeki NATP00295/NATP00296 kesitlerinden
  doğrulanamadı. NATP00296 Problem 2, satır 184'teki ters işlemi doğrudan işlemle denetleme pasajı doğrulandı.

## Teslim sırası

1. `71438dd`: izole Hermite/konik modülü ve ilk birim testleri; main'e gönderildi.
2. `c645594`: dil entegrasyonu, sayısal değerler, parametre ve kimlik düzeltmeleri; main'e gönderildi.
3. `1d7c06a`: test fixture'ı ve bc regresyonları; main'e gönderildi.
4. README, karar kaydı, handoff ve bu rapor ayrı belge commitinde yer alıyor.
