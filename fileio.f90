subroutine read_restart_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
                             linkage_block_fitness,           &
                             initial_allele_effects,          &
                             generation_number,max_size,myid_str) 

! This routine reads a dump file containing the mutation arrays
! dmutn and fmutn and the linkage block mutation count array
! lb_mutn_count and the linkage block fitness array 
! linkage_block_fitness for purposes of restart.  The argument
! generation_number is the generation number of the dump file 
! being read.
use inputs
use profiler
include 'common.h'
integer generation_number, max_size, i, lb, n1, n2, npath, nm
integer dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer nmutn(max_neu_mutn_per_indiv/2,2,max_size)
integer fmutn(max_fav_mutn_per_indiv/2,2,max_size)
integer lb_mutn_count(num_linkage_subunits,2,3,max_size)
real*8 linkage_block_fitness(num_linkage_subunits,2,max_size)
real  initial_allele_effects(num_linkage_subunits)
real f1
logical l1
character char1*1, char12*12, char20*20, char32*32, path*80
character myid_str*3
call second(tin)

npath = index(data_file_path,' ') - 1
char1 = char(48 + restart_dump_number)

open (10, file=data_file_path(1:npath)//case_id// &
      '.'//myid_str//'.dmp.'//char1,status='unknown')

read(10,'(a20,i12)') char32, generation_number

dmutn = num_linkage_subunits*lb_modulo + 1
fmutn = num_linkage_subunits*lb_modulo + 1

do i=1,pop_size
      read(10,*) ! read label '=== individual: # ==='
      read(10,*) ! read label 'lb_mutn_count:'
      read(10,'(12i6)') lb_mutn_count(:,:,:,i) 
      read(10,*) ! read label 'linkage_block_fitness:'
      read(10,'(6f12.8)') linkage_block_fitness(:,:,i) 

      read(10,*) ! read label 'deleterious mutations:'
      read(10,'( i12)') dmutn(1,1,i)
      nm = min(max_del_mutn_per_indiv/2, dmutn(1,1,i)+1)
      read(10,'(6i12)') dmutn(2:nm,1,i)
      read(10,'( i12)') dmutn(1,2,i)
      nm = min(max_del_mutn_per_indiv/2, dmutn(1,2,i)+1)
      read(10,'(6i12)') dmutn(2:nm,2,i)

      read(10,*) ! read label 'neutral mutations:'
      read(10,'( i12)') nmutn(1,1,i)
      nm = min(max_neu_mutn_per_indiv/2, nmutn(1,1,i)+1)
      read(10,'(6i12)') nmutn(2:nm,1,i)
      read(10,'( i12)') nmutn(1,2,i)
      nm = min(max_neu_mutn_per_indiv/2, nmutn(1,2,i)+1)
      read(10,'(6i12)') nmutn(2:nm,2,i)

      read(10,*) ! read label 'favorable mutations:'
      read(10,'( i12)') fmutn(1,1,i)
      nm = min(max_fav_mutn_per_indiv/2, fmutn(1,1,i)+1)
      read(10,'(6i12)') fmutn(2:nm,1,i)
      read(10,'( i12)') fmutn(1,2,i)
      nm = min(max_fav_mutn_per_indiv/2, fmutn(1,2,i)+1)
      read(10,'(6i12)') fmutn(2:nm,2,i)
end do

if(num_contrasting_alleles > 0) then
   read(10,'(6f12.9)') initial_allele_effects
   if(.not. is_parallel) &
   write(6, '("Restart run will use the previous value for"/ &
              "parameter num_contrasting_alleles =", i10)')  &
                         num_contrasting_alleles
   write(9, '("Restart run will use the previous value for"/ &
              "parameter num_contrasting_alleles =", i10)')  &
                         num_contrasting_alleles
end if

close (10)

call second(tout)
sec(9) = sec(9) + tout - tin
end subroutine read_restart_dump

subroutine read_mutn_file(dmutn,nmutn,fmutn,lb_mutn_count, &
                          linkage_block_fitness,max_size)
! This routine reads the file caseid_mutn.in only when upload_mutations
! flag is set to 1.
use inputs
include 'common.h'
integer id, lb, hap_id, mutn, dominance
real    fitness, w
integer i, npath, max_size
integer nimpi(pop_size), encode_mutn
integer dmutn(max_del_mutn_per_indiv/2,2,max_size)
integer nmutn(max_neu_mutn_per_indiv/2,2,max_size)
integer fmutn(max_fav_mutn_per_indiv/2,2,max_size)
integer lb_mutn_count(num_linkage_subunits,2,3,max_size)
real*8 linkage_block_fitness(num_linkage_subunits,2,max_size)

npath = index(data_file_path,' ') - 1

open (10, file=data_file_path(1:npath)//case_id// &
      '_mutn.in',status='unknown')

read(10,*) num_uploaded_mutn
write(*,*) 'Reading mutation file with ', num_uploaded_mutn, &
           'mutations'
write(*,*) '#     id, linkage_block, haplotype, fitness, dominant/recessive, mutn_id'
read(10,*) ! header

w = multiplicative_weighting

do i=1,num_uploaded_mutn

   read(10,*) id,lb,hap_id,fitness,dominance
   mutn = encode_mutn(fitness,lb,dominance)
   uploaded_mutn(i) = mutn
   write(*,*) id,lb,hap_id,fitness,dominance,mutn

   if(fitness > 0.) then
      fmutn(1,hap_id,id) = fmutn(1,hap_id,id) + 1
      fmutn(fmutn(1,hap_id,id)+1,hap_id,id) = mutn
      lb_mutn_count(lb,hap_id,2,id) = lb_mutn_count(lb,hap_id,2,id) + 1
   elseif(fitness < 0.) then
      dmutn(1,hap_id,id) = dmutn(1,hap_id,id) + 1
      dmutn(dmutn(1,hap_id,id)+1,hap_id,id) = mutn
      lb_mutn_count(lb,hap_id,1,id) = lb_mutn_count(lb,hap_id,1,id) + 1
   else ! fitness = 0 --> neutral mutation
      nmutn(1,hap_id,id) = nmutn(1,hap_id,id) + 1
      nmutn(nmutn(1,hap_id,id)+1,hap_id,id) = mutn
      lb_mutn_count(lb,hap_id,3,id) = lb_mutn_count(lb,hap_id,3,id) + 1
   end if

   linkage_block_fitness(lb,hap_id,id) =   &
      (linkage_block_fitness(lb,hap_id,id) &
       + (1. - w)*fitness) * (1.d0 + w*fitness)
  
end do
   
close (10)

end subroutine read_mutn_file
     
subroutine write_output_dump(dmutn,nmutn,fmutn,lb_mutn_count, &
                             linkage_block_fitness,     &
                             initial_allele_effects,    &
                             generation_number,myid_str)

! This routine writes an output file containing the current
! parameter values, the stored mutation arrays dmutn and fmutn,
! the linkage block mutation count array lb_mutn_count and
! the linkage block fitness array linkage_block_fitness for
! the current generation specified by generation_number.
use inputs
use profiler
include 'common.h'
integer dmutn(max_del_mutn_per_indiv/2,2,*), &
        fmutn(max_fav_mutn_per_indiv/2,2,*), &
        nmutn(max_neu_mutn_per_indiv/2,2,*), &
               lb_mutn_count(num_linkage_subunits,2,3,*) 
real*8 linkage_block_fitness(num_linkage_subunits,2,*)
real  initial_allele_effects(num_linkage_subunits)
integer generation_number, i, lb, npath
character char1*1
character myid_str*3
call second(tin)

npath = index(data_file_path,' ') - 1
char1 = char(48 + dump_number)

open (10, file=data_file_path(1:npath)//case_id// &
      '.'//myid_str//'.dmp.'//char1,status='unknown')

write(10,'(a20, i12)') 'generation_number = ', generation_number

do i=1,pop_size
      write(10,*) '=========== individual:',i,'================================='

      write(10,*) 'lb_mutn_count:'
      write(10,'(12i6)') lb_mutn_count(:,:,:,i) 
      write(10,*) 'linkage_block_fitness:'
      write(10,'(6f12.8)') linkage_block_fitness(:,:,i) 

      write(10,*) 'deleterious mutations:'
      write(10,'( i12)') dmutn(1,1,i)
      write(10,'(6i12)') dmutn(2:dmutn(1,1,i)+1,1,i)
      write(10,'( i12)') dmutn(1,2,i)
      write(10,'(6i12)') dmutn(2:dmutn(1,2,i)+1,2,i)

      write(10,*) 'neutral mutations:'
      write(10,'( i12)') nmutn(1,1,i)
      write(10,'(6i12)') nmutn(2:nmutn(1,1,i)+1,1,i)
      write(10,'( i12)') nmutn(1,2,i)
      write(10,'(6i12)') nmutn(2:nmutn(1,2,i)+1,2,i)

      write(10,*) 'favorable mutations:'
      write(10,'( i12)') fmutn(1,1,i)
      write(10,'(6i12)') fmutn(2:fmutn(1,1,i)+1,1,i)
      write(10,'( i12)') fmutn(1,2,i)
      write(10,'(6i12)') fmutn(2:fmutn(1,2,i)+1,2,i)
end do

if(num_contrasting_alleles > 0)  &
   write(10,'(6f12.9)') initial_allele_effects

close (10)

call second(tout)
sec(10) = sec(10) + tout - tin
end subroutine write_output_dump

subroutine write_vcf_file(dmutn)

use inputs
use polygenic
include 'common.h'

integer :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer :: npath, i, j, k, lb, dominance
integer :: chrom, lb_per_chrom, pos
real*8 :: x, fitness, id
character*8 :: pos_str
character*7 :: id_str
character*5 :: chrom_str
character*3 :: myid_str
character*1 :: dot = ".", ref, alt
character*1 :: tab = char(9)
integer, external :: decode_mutn_id

npath = index(data_file_path,' ') - 1

if (is_parallel) then
   write(myid_str,'(i3.3)') myid+1
else
   write(myid_str,'(i3.3)') 0
   myid = 0
end if

open(27, file=data_file_path(1:npath)//case_id//'.'//myid_str &
//'.vcf',status='unknown')

write(27,'("##fileformat=VCFv4.1")')
write(27,'("##contig=<ID=1,length=249250621,assembly=b37>")')
write(27,'(a, a)') "####reference=", case_id
write(27,'("#CHROM",a1,"POS",a1,"ID",a1,"REF",a1,"ALT",a1,"QUAL",a1,"FILTER",a1,"INFO")') &
           tab, tab, tab, tab, tab, tab, tab

lb_per_chrom = num_linkage_subunits/haploid_chromosome_number

do k = 1, pop_size
   do j = 1, 2
      do i = 2, dmutn(1,j,k) 
         call decode_mutn_del(dmutn(i,j,k), lb, dominance, fitness)
         if (lb == 0) cycle ! if for some reason the linkage block number
                            ! is zero ignore for now... for some reason about
                            ! 3 out of every 11,000 show up as zero... not sure why
         x = lb/real(lb_per_chrom)
         chrom = ceiling(x)
         pos = ceiling((x - chrom + 1)*lb_per_chrom)*lb_modulo
         !id = decode_mutn_id(dmutn(i,j,k), -1) ! don't know why I'm having so much trouble w/this!!
         id = real(mod(dmutn(i,j,k), lb_modulo))*del_scale
         ref = "X" !random_nucl()
         alt = mutn_to_nucl(dmutn(i,j,k))
         do while(ref == alt) ! generate a value of alt that is different than ref
            alt = random_nucl()
         end do
         write(chrom_str, '(a, i2.2)') "chr", chrom
         write(pos_str, '(i0.8)') pos
         write(id_str, '(f7.5)') id
         write(27,'(15a)') chrom_str//tab//pos_str//tab//id_str//tab//ref//tab//alt//tab//dot//tab//dot//tab//dot
      end do
   end do 
end do

close(27)

end subroutine write_vcf_file

subroutine write_alleles(dmutn)

use inputs
include 'common.h'
integer :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer :: id, i, j, k, npath
integer, parameter :: intmax = 2147483647
character*3 :: myid_str
character*1 :: comma = ","

npath = index(data_file_path,' ') - 1

if (is_parallel) then
   write(myid_str,'(i3.3)') myid+1
else
   write(myid_str,'(i3.3)') 0
   myid = 0
end if

!open(27, file=data_file_path(1:npath)//case_id//'.'//myid_str &
!//'.csv',status='unknown')
open(27, file='alleles.csv', status='unknown')

do k = 1, pop_size
   write(27, '(i12, a1, $)') k, comma
   do j = 1, 2
      do i = 2, dmutn(1,j,k) 
         id = int(intmax * real(mod(dmutn(i,j,k), lb_modulo))*del_scale)
         write(27, '(i12, a1, $)') dmutn(i,j,k), comma
      end do
   end do
   write(27, *)
end do

close(27)

end subroutine write_alleles

subroutine write_sample(dmutn,fmutn,lb_mutn_count, &
                        linkage_block_fitness,fitness, &
                        defect,improve,effect,del_mutn, &
                        fav_mutn,generation_number)

! This routine writes an output file containing details concerning
! the mutations carried by five members of the total population in
! a format that is (hopefully) readily understandable.
use inputs
use profiler
include 'common.h'
integer :: dmutn(max_del_mutn_per_indiv/2,2,*)
integer :: fmutn(max_fav_mutn_per_indiv/2,2,*)
integer :: lb_mutn_count(num_linkage_subunits,2,2,*)
integer :: del_mutn(num_linkage_subunits,2)
integer :: fav_mutn(num_linkage_subunits,2)
real*8  :: linkage_block_fitness(num_linkage_subunits,2,*)
real*8  :: fitness(*), x
real    :: defect(num_linkage_subunits,2,*)
real    :: improve(num_linkage_subunits,2,*)
real    :: effect(*)
integer :: generation_number
integer, parameter :: fid = 20

integer :: i, j, lb, m, mutn
real    :: d
call second(tin)

rewind (fid)

call write_parameters(fid)

write(fid,'(/23x,"generation number = ",i6)') generation_number

do i=current_pop_size/10,current_pop_size/3,current_pop_size/5

   write(fid,'(/"individual number = ",i6,"  fitness = ",f9.6)') &
         i, fitness(i)

   del_mutn = 0
   fav_mutn = 0

   do m=2,dmutn(1,1,i)+1
      lb = abs(dmutn(m,1,i))/lb_modulo + 1
      del_mutn(lb,1) = del_mutn(lb,1) + 1
      x  = mod(abs(dmutn(m,1,i)), lb_modulo)*del_scale
      d  = dexp(-alpha_del*x**gamma_del)
      if(x >= 1.d0) d = 0.
      if(dmutn(m,1,i) < 0) d = -d
      defect(lb,1,del_mutn(lb,1)) = d
   end do

   do m=2,dmutn(1,2,i)+1
      lb = abs(dmutn(m,2,i))/lb_modulo + 1
      del_mutn(lb,2) = del_mutn(lb,2) + 1
      x  = mod(abs(dmutn(m,2,i)), lb_modulo)*del_scale
      d  = dexp(-alpha_del*x**gamma_del)
      if(x >= 1.d0) d = 0.
      if(dmutn(m,2,i) < 0) d = -d
      defect(lb,2,del_mutn(lb,2)) = d
   end do

   do m=2,fmutn(1,1,i)+1
      lb = abs(fmutn(m,1,i))/lb_modulo + 1
      fav_mutn(lb,1) = fav_mutn(lb,1) + 1
      mutn = mod(abs(fmutn(m,1,i)), lb_modulo)
      d    = dexp(-alpha_fav*(real(mutn)*fav_scale)**gamma_fav) &
             *max_fav_fitness_gain
      if(fmutn(m,1,i) < 0) d = -d
      improve(lb,1,fav_mutn(lb,1)) = d
   end do

   do m=2,fmutn(1,2,i)+1
      lb = abs(fmutn(m,2,i))/lb_modulo + 1
      fav_mutn(lb,2) = fav_mutn(lb,2) + 1
      mutn = mod(abs(fmutn(m,2,i)), lb_modulo)
      d    = dexp(-alpha_fav*(real(mutn)*fav_scale)**gamma_fav) &
             *max_fav_fitness_gain
      if(fmutn(m,2,i) < 0) d = -d
      improve(lb,2,fav_mutn(lb,2)) = d
   end do

   do lb=1,num_linkage_subunits

      write(fid,'(/27x,"lb number = ",i6)') lb

      write(fid,'("Haplotype 1: total deleterious mutn count = ", &
                 i6, "  composite fitness = ",f9.6)')  &
         lb_mutn_count(lb,1,1,i), linkage_block_fitness(lb,1,i)

      if(tracking_threshold /= 1.0) then
  
      j = 0
      do m=1,del_mutn(lb,1)
         if(defect(lb,1,m) < 0) then
            j = j + 1
            effect(j) = -defect(lb,1,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness degradations of tracked deleterious " &
                    "recessive mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,del_mutn(lb,1)
         if(defect(lb,1,m) > 0) then
            j = j + 1
            effect(j) = defect(lb,1,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness degradations of tracked deleterious " &
                    "dominant mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,fav_mutn(lb,1)
         if(improve(lb,1,m) < 0) then
            j = j + 1
            effect(j) = -improve(lb,1,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness improvements of tracked favorable " &
                    "recessive mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,fav_mutn(lb,1)
         if(improve(lb,1,m) > 0) then
            j = j + 1
            effect(j) = improve(lb,1,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness improvements of tracked favorable " &
                    "dominant mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      end if 

      write(fid,'("Haplotype 2: total deleterious mutn count = ", &
                 i6, "  composite fitness = ",f9.6)') &
         lb_mutn_count(lb,2,1,i), linkage_block_fitness(lb,2,i) 

      if(tracking_threshold /= 1.0) then

      j = 0
      do m=1,del_mutn(lb,2)
         if(defect(lb,2,m) < 0) then
            j = j + 1
            effect(j) = -defect(lb,2,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness degradations of tracked deleterious " &
                    "recessive mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,del_mutn(lb,2)
         if(defect(lb,2,m) > 0) then
            j = j + 1
            effect(j) = defect(lb,2,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness degradations of tracked deleterious " &
                    "dominant mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,fav_mutn(lb,2)
         if(improve(lb,2,m) < 0) then
            j = j + 1
            effect(j) = -improve(lb,2,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness improvements of tracked favorable " &
                    "recessive mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      j = 0
      do m=1,fav_mutn(lb,2)
         if(improve(lb,2,m) > 0) then
            j = j + 1
            effect(j) = improve(lb,2,m)
         end if
      end do

      if(j > 0) then
         write(fid,'("Fitness improvements of tracked favorable " &
                    "dominant mutations:")')
         write(fid,'(8f9.6)') (effect(m),m=1,j)
      end if 

      end if 

   end do

end do

call flush(fid)

call second(tout)
sec(11) = sec(11) + tout - tin

end subroutine write_sample

subroutine write_status(unit, gen, current_pop_size, &
           frac_recessive, total_del_mutn, tracked_del_mutn, &
           total_fav_mutn, neu_mutn, pre_sel_fitness, &
           pre_sel_geno_sd, pre_sel_pheno_sd, pre_sel_corr, &
           post_sel_fitness, post_sel_geno_sd, post_sel_pheno_sd, &
           post_sel_corr, num_polys_this_gen, num_polys_cumulative )
use inputs
integer unit, gen, current_pop_size, num_polygenics
real*8  pre_sel_fitness, pre_sel_geno_sd, pre_sel_pheno_sd, &
        pre_sel_corr, post_sel_fitness, post_sel_geno_sd, &
        post_sel_pheno_sd, post_sel_corr 
real*8  total_del_mutn, tracked_del_mutn, total_fav_mutn, &
        frac_recessive, neu_mutn

!call system('clear') ! clear the screen before each output
write(unit,'(/"generation =",i10,"  population size =", i6, &
      "  frac recessive =",f7.4/"before sel: geno fitness =",f9.5, &
      "  geno s.d. =",f8.5,"  pheno s.d. =",f8.5/ &
      "after  sel:               ",f9.5,"             ",f8.5, &
      "              ",f8.5/"before sel geno-pheno corr =",f7.4, &
      "    after sel geno-pheno corr  =",f10.4/ &
      "del mutn/indiv =",i6, &
      "  tracked del/ind =",i6,"  fav mutn/indiv   =",f10.4)') &
      gen, current_pop_size, frac_recessive, &
      pre_sel_fitness, pre_sel_geno_sd, pre_sel_pheno_sd, &
      post_sel_fitness, post_sel_geno_sd, post_sel_pheno_sd, &
      pre_sel_corr, post_sel_corr, &
      int(total_del_mutn/current_pop_size), &
      int(tracked_del_mutn/current_pop_size), &
      total_fav_mutn/current_pop_size 
      if(neu_mutn > 0 .and. .not.polygenic_beneficials) write(unit, &
         '("neu mutn/indiv =",i6,$)'), int(neu_mutn/current_pop_size)
      if(num_polys_cumulative > 0) then
         write(unit, '("polys_this_gen =",i6,2x,"polys_cumulativ =",i6)') &
               num_polys_this_gen, num_polys_cumulative
      else
         write(unit, *)
      endif 

end subroutine write_status
