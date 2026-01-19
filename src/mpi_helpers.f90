module mpi_helpers
use mpi

interface mpi_davg
   module procedure mpi_davg_r1
   module procedure mpi_davg_r2
end interface

interface mpi_dsum
   module procedure mpi_dsum_r0
   module procedure mpi_dsum_r1
   module procedure mpi_dsum_r2
end interface

interface mpi_ravg
   module procedure mpi_ravg_r0
   module procedure mpi_ravg_r1
   module procedure mpi_ravg_r2
end interface

interface mpi_isum
   module procedure mpi_isum_r0
   module procedure mpi_isum_r1
   module procedure mpi_isum_r2
end interface

interface mpi_mybcastd
   module procedure mpi_mybcastd_r0
   module procedure mpi_mybcastd_r1
   module procedure mpi_mybcastd_r2
end interface

interface mpi_mybcasti
   module procedure mpi_mybcasti_r0
   module procedure mpi_mybcasti_r1
   module procedure mpi_mybcasti_r2
end interface

interface mpi_gather_d
   module procedure mpi_gather_d_r0
   module procedure mpi_gather_d_r1
end interface

contains

subroutine mpi_myinit(ownid,ierr_out)
   include 'common.h'
   integer ownid, ierr_out
   call MPI_INIT(ierr_out)
   call MPI_COMM_DUP(MPI_COMM_WORLD,MYCOMM,ierr_out)
   call MPI_COMM_RANK(MYCOMM,ownid,ierr_out)
   call MPI_COMM_SIZE(MYCOMM,num_tribes,ierr_out)
end

subroutine mpi_mybarrier()
   include 'common.h'
   call MPI_BARRIER(MYCOMM,ierr)
end

subroutine mpi_myabort()
   include 'common.h'
   call MPI_Abort(MYCOMM,1,ierr)
end

subroutine mpi_mybcastd_r0(x,n)
   include 'common.h'
   integer n
   real*8 x
   call MPI_Bcast(x,n,MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)
end

subroutine mpi_mybcastd_r1(x,n)
   include 'common.h'
   integer n
   real*8 x(:)
   call MPI_Bcast(x,n,MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)
end

subroutine mpi_mybcastd_r2(x,n)
   include 'common.h'
   integer n
   real*8 x(:,:)
   call MPI_Bcast(x,n,MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)
end

subroutine mpi_mybcasti_r0(x,n)
   include 'common.h'
   integer x, n
   call MPI_Bcast(x,n,MPI_INTEGER,0,MYCOMM,ierr)
end

subroutine mpi_mybcasti_r1(x,n)
   include 'common.h'
   integer n
   integer x(:)
   call MPI_Bcast(x,n,MPI_INTEGER,0,MYCOMM,ierr)
end

subroutine mpi_mybcasti_r2(x,n)
   include 'common.h'
   integer n
   integer x(:,:)
   call MPI_Bcast(x,n,MPI_INTEGER,0,MYCOMM,ierr)
end

subroutine mpi_send_int(istatus,dest,msg_num_in,ierr_out)
   include 'common.h'
   integer istatus,dest,msg_num_in,ierr_out
   call MPI_Send(istatus,1,MPI_INTEGER,dest,msg_num_in,MYCOMM,ierr_out)
end

subroutine mpi_recv_int(istatus,src,msg_num_in,ierr_out)
   include 'common.h'
   integer status(MPI_STATUS_SIZE)
   integer istatus,src,msg_num_in,ierr_out
   call MPI_Recv(istatus,1,MPI_INTEGER,src,msg_num_in,MYCOMM,status,ierr_out)
end

subroutine mpi_myfinalize(ierr_out)
   integer ierr_out
   call MPI_FINALIZE(ierr_out)
end

subroutine mpi_davg_scalar(x, xavg, n)
   include 'common.h'
   integer n
   real*8 x, xavg
   call MPI_REDUCE(x,xavg,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_davg_r1(x, xavg, n)
   include 'common.h'
   integer n
   real*8 x(:), xavg(:)
   call MPI_REDUCE(x,xavg,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_davg_r2(x, xavg, n)
   include 'common.h'
   integer n
   real*8 x(:,:), xavg(:,:)
   call MPI_REDUCE(x,xavg,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_dsum_r0(x, xsum, n)
   include 'common.h'
   integer n
   real*8 x, xsum
   call MPI_REDUCE(x,xsum,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_dsum_r1(x, xsum, n)
   include 'common.h'
   integer n
   real*8 x(:), xsum(:)
   call MPI_REDUCE(x,xsum,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_dsum_r2(x, xsum, n)
   include 'common.h'
   integer n
   real*8 x(:,:), xsum(:,:)
   call MPI_REDUCE(x,xsum,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_ravg_r0(x, xavg, n)
   include 'common.h'
   integer n
   real x, xavg
   call MPI_REDUCE(x,xavg,n,MPI_REAL,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_ravg_r1(x, xavg, n)
   include 'common.h'
   integer n
   real x(:), xavg(:)
   call MPI_REDUCE(x,xavg,n,MPI_REAL,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_ravg_r2(x, xavg, n)
   include 'common.h'
   integer n
   real x(:,:), xavg(:,:)
   call MPI_REDUCE(x,xavg,n,MPI_REAL,MPI_SUM,0,MYCOMM,ierr)
   xavg = xavg/num_tribes
end

subroutine mpi_isum_r0(i, isum, n)
   include 'common.h'
   integer i, isum, n
   call MPI_REDUCE(i,isum,n,MPI_INTEGER,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_isum_r1(i, isum, n)
   include 'common.h'
   integer i(:), isum(:), n
   call MPI_REDUCE(i,isum,n,MPI_INTEGER,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_isum_r2(i, isum, n)
   include 'common.h'
   integer i(:,:), isum(:,:), n
   call MPI_REDUCE(i,isum,n,MPI_INTEGER,MPI_SUM,0,MYCOMM,ierr)
end

subroutine mpi_gather_d_r0(sendbuf, recvbuf, n)
   include 'common.h'
   integer n
   real*8 sendbuf
   real*8 recvbuf(:)
   call MPI_GATHER(sendbuf,n,MPI_DOUBLE_PRECISION, &
                   recvbuf,n,MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)
end

subroutine mpi_gather_d_r1(sendbuf, recvbuf, n)
   include 'common.h'
   integer n
   real*8 sendbuf(:)
   real*8 recvbuf(:)
   call MPI_GATHER(sendbuf,n,MPI_DOUBLE_PRECISION, &
                   recvbuf,n,MPI_DOUBLE_PRECISION,0,MYCOMM,ierr)
end

end module mpi_helpers
