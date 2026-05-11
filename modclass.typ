#import "@preview/mitex:0.3.1": mi, mitex

#heading(level: 1)[สรุปการดัดแปลงโค้ด CLASS สำหรับแบบจำลอง Quintom]

การดัดแปลงโค้ด Boltzmann (CLASS) เพื่อรองรับแบบจำลองพลังงานมืดแบบ Quintom (ผสมผสานระหว่าง Quintessence และ Phantom field) รวมถึงอันตรกิริยากับสสาร (Coupling) ประกอบไปด้วยการแก้ไขโครงสร้างโค้ดหลักๆ ในภาษา C และส่วนเชื่อมต่อในภาษา Python ดังรายละเอียดต่อไปนี้:

#heading(level: 2)[1. โครงสร้างตัวแปรของสมการพื้นหลัง (Background Variables)]
มีการเพิ่มตัวแปรใหม่เข้าไปใน `struct background` ของไฟล์ `include/background.h` และส่วน C-interface (`python/cclassy.pxd`) ให้รองรับพารามิเตอร์ของแบบจำลอง Quintom (`qtm`) อย่างครบถ้วน ได้แก่:
- *พารามิเตอร์ความหนาแน่นตั้งต้น:* `Omega0_qtm` (ความหนาแน่นสัมพัทธ์), `rho0_qtm`, `p0_qtm`
- *พารามิเตอร์ของศักย์ (Potential) และอันตรกิริยา (Coupling):* `lambda_qtm` (ความชันของศักย์ #mi(`\lambda`)), `delta_qtm` (ค่าคงที่อันตรกิริยาระหว่าง Phantom และสสาร #mi(`\delta`)), `logV0_qtm` (แอมพลิจูดของศักย์ #mi(`\ln V_0`))
- *พารามิเตอร์สเกลของสสารคู่ควบ:* `Kb_qtm` (สำหรับบาริออน) และ `Kcdm_qtm` (สำหรับสสารมืดเย็น)
- *แฟลกการตั้งค่า:* `has_qtm` ตรวจสอบการเปิดใช้ Quintom รวมถึงแฟลก `coupled_baryon_qtm` และ `coupled_cdm_qtm` เพื่อระบุว่าจะให้อนุภาคใดมีอันตรกิริยากับสนาม Phantom
- *ดัชนี (Indices):* นำดัชนีใหม่เข้าสู่สมการสถานะพื้นหลังสำหรับการอินทิเกรต ได้แก่ สนาม Quintessence #mi(`\phi, \phi'`), สนาม Phantom #mi(`\sigma, \sigma'`), สมการศักย์ #mi(`V, \partial V, \partial^2 V`) และความหนาแน่นส่วนพลังงานกับความดัน

#heading(level: 2)[2. สมการวิวัฒนาการและฟิสิกส์พื้นหลัง (Background Physics)]
การเปลี่ยนแปลงสำคัญที่สุดในฝั่งสมการเชิงอนุพันธ์เกิดขึ้นในไฟล์ `source/background.c`:
- *ฟังก์ชันศักย์:* นิยามฟังก์ชัน `V_qtm`, `dV_qtm`, และ `ddV_qtm` เพื่อสร้างศักย์ #mi(`V(\phi) = V_0 e^{-\lambda \phi}`) โดยผูกกับสนาม Quintessence
- *เงื่อนไขเริ่มต้น (Initial Conditions):* ถูกตั้งค่าใน `background_initial_conditions` ให้สนาม #mi(`\phi`) เริ่มที่สภาวะแกะรอย (Attractor solution) ส่วนสนาม #mi(`\sigma`) และ #mi(`\sigma'`) เริ่มต้นที่ 0 
- *สมการวิวัฒนาการความหนาแน่น (Density & Pressure):* ใน `background_functions` ได้กำหนดให้ส่วน Quintom มีทั้งแบบ Quintessence พลังงานจลน์บวก และ Phantom พลังงานจลน์ลบ (Negative kinetic energy):
  #mitex(`\rho_\sigma = - \frac{(\sigma')^2}{6 a^2}, \quad p_\sigma = - \frac{(\sigma')^2}{6 a^2}`)
  ขณะเดียวกัน หากมีการเปิดใช้ Coupling ค่าความหนาแน่นของสสาร (CDM และ/หรือ แบริออน) จะเปลี่ยนรูปแบบไปขึ้นกับการคู่ควบ:
  #mitex(`\rho_m = K_{m} e^{\delta \sigma} a^{-3}`)
- *สมการการเคลื่อนที่ (Klein-Gordon equations):* จัดการอนุพันธ์ย่อยใน `background_derivs`:
  - สนาม #mi(`\phi`): 
    #mitex(`\phi'' + 2aH \phi' = -a^2 \frac{\partial V}{\partial \phi}`)
  - สนาม #mi(`\sigma`): มีพจน์ Source term จากสสารมาร่วมด้วยตามค่าทฤษฎี: 
    #mitex(`\sigma'' + 2aH \sigma' = 3a \delta \rho_{M} / H`)

#heading(level: 2)[3. ระเบียบวิธี Shooting (Input & Shooting Parameters)]
เพื่อที่จะให้ความหนาแน่นรวมสอดคล้องกับเป้าหมายที่ตั้งไว้ในจักรวาลวิทยา โค้ดใน `source/input.c` ถูกแก้ไขให้ยิงค่า (Shooting Method) ได้:
- *เป้าหมายการยิง:* เพิ่ม `Omega_qtm`, `omega_b_qtm`, `omega_cdm_qtm` ให้เลือกใช้เป้าหมายในการแก้สมการ 
- *ค่าตั้งต้น (Guess):* สำหรับ `Omega_qtm` ตัวยิงถูกจับคู่กับค่าตั้งต้น `logV0_qtm` (โดยกะเส้นแนวโน้มจาก #mi(`\ln(\Omega_{qtm} H_0^2)`)) ในขณะที่ของสสารตั้งเป็นตัวแปรสเกล #mi(`K_b`) และ #mi(`K_{cdm}`) 
- มีการรวมค่าเป้าหมายของ `Omega_qtm` เข้ากับกฎข้อขัดแย้งเชิงความหนาแน่น (Budget equation constraints) ฝั่งพารามิเตอร์อื่นๆ ด้วย

#heading(level: 2)[4. ส่วนเชื่อมต่อภาษาไพธอน (Python Wrapper)]
เพื่อให้นำมาวิเคราะห์ในรูป MCMC หรือผ่าน Jupyter Notebook ได้สะดวก มีการขยายคลาสใน `python/classy.pyx`:
- เพิ่ม properties อย่าง `@property def Omega_qtm`, `lambda_qtm`, `delta_qtm`, `w_qtm`, ฯลฯ ลงในโครงสร้างคลาส `Class` ของ Python
- ค่าสมการสถานะ (EoS) ของสสารมืด Quintom #mi(`w_{qtm}`) ในปัจจุบันได้รับการคำนวณผ่าน #mi(`p_0 / \rho_0`) ส่งออกมาให้บริการแก่นักพัฒนาและสคริปต์ทางสถิติ
  
สรุปได้ว่าการปรับเปลี่ยนนี้ขยายขีดความสามารถของรหัส CLASS เดิม ให้ครอบคลุมการทำงานทั้งการวิวัฒนาการของสมการอนุพันธ์แบบจำลอง Quintom และการส่งอินพุตรับพารามิเตอร์อย่างเต็มรูปแบบ
