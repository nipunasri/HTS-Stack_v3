#define _OCL_CODE_
#include "HTSShared.hpp"

#pragma OPENCL EXTENSION cl_amd_printf:enable // enable printf on AMD GPU's

/************************* uiAlloc *****************************/
/* Function name	: uiAlloc
 * Arguments		: NodePool, pointer-to-free-node in NodePool	
 * Return value		: Index of a free node in nodepool
 * Description		: This function returns index of a free node, whose freebit is 1 (3rd bit), 
					  A work-item level function.
 */
int uiAlloc( TLLNode *pNodePool, /* Hash-table + stack of free nodes */ 
			  uint uiHeadIndex ) /* pointer to free-node in pNodePool */
{ 
	// private variables 
    uint uiLLMNode = 0 ; /* [<next-node-index>,<Bits(f,r,m)>] */
    int uiLLNode = 0 ; /* next-node-index in nodepool, shift 3 bits right to get */ 
    bool bFoundNode =  false  ; /* true, node found otherwise not found */
    uint freeBit = 0 ;

    while( bFoundNode != true ) { /* iterate untill NodeFound */

		uiLLMNode = pNodePool[uiHeadIndex].uiNext ; /* Always top node is free in stack */
		uiLLNode = GET_PTR(uiLLMNode) ;/* Get actual index of next node in pool */
            
		if(uiLLNode >= OCL_NODE_POOL_SIZE) return NULL ; // BC

		 uint uiLLMNextNode = pNodePool[uiLLNode].uiNext ; /* Get next node index */
		 uint uiLLNextNode = GET_PTR(uiLLMNextNode) ;
			freeBit = GET_FBIT(uiLLMNextNode) ; /* Get free Bit */
			if(freeBit) { /* If, this is a free node */	
            atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[uiHeadIndex].uiNext));
            // update starting index
            bFoundNode = atomic_compare_exchange_strong
                (pChgPtr,
                &uiLLMNode,
                uiLLMNextNode); // POINTS TO NEXT FREE NODE
            } // end-of-if
     } // end-of-while	
    // clear the data, if any
    for(int i = 0; i < OCL_WG_SIZE; i++) 
		pNodePool[uiLLNode].pE[i] = 0 ;
	pNodePool[uiLLNode].uiNext = 0 ;
	return uiLLNode ; 
} // end-of-uiAlloc

/************************* uiFree *****************************/
/* Function name	: uiFree
 * Arguments		: NodePool(pointer to nodepool), Node-index to be deleted	
 * Return value		: void
 * Description		: This function marks given index node as free node, by setting free bit,
					  executes at work-item level
 */
bool uiFree( TLLNode *pNodePool,	/* Hash-table + stack of free nodes */
			  uint uiHeadIndex,		/* pointer to free-node in pNodePool */
			  uint iDelNode )		/* node to be deleted */	
{ 
    uint uiLLMNode = 0 ; /* [<next-node-index>,<Bits(f,r,m)>] */
    uint uiLLNode = 0 ; /* next-node-index in nodepool, shift 3 bits right to get */ 
    bool bFoundNode =  false  ; /* true, node found otherwise not found */
    uint uiLLDelMNode = 0 ;
    uint uiLLNewNode = 0 ;
    uint freeBit = 0 ;
    
    uiLLDelMNode = pNodePool[iDelNode].uiNext ;
    freeBit = GET_FBIT(uiLLDelMNode) ;
    if(freeBit) return bFoundNode ;

    while( bFoundNode != true ) { /* iterate untill NodeFound */

        uiLLMNode = pNodePool[uiHeadIndex].uiNext ; /* Always top node is free : stack */
        uiLLNode = GET_PTR(uiLLMNode) ; /* Get actual index of next node in pool */
        uiLLDelMNode = SET_FBIT(uiLLMNode); /* index with free bit set */

        pNodePool[iDelNode].uiNext = uiLLDelMNode; // points to first node

        atomic_uint* pChgPtr =
            (atomic_uint *)(&(pNodePool[uiHeadIndex].uiNext));

        uiLLNewNode = SET_PTR(iDelNode) ;
        
        bFoundNode = atomic_compare_exchange_strong
            (pChgPtr,
            &uiLLMNode,
            uiLLNewNode);
        } // end-of-while-loop
    return bFoundNode ;
} // end-of-free

/************************* SNIP *****************************/
/* Function name : SNIP
 * Arguments	 : pNodePool, HeadIndex
 * Description   : This function takes starting index as input and deletes next marked node,
                   if it is marked for deletion otherwise return next-node-un-marked.
 * Return values : Returns snipped on successfull deletion of marked node, otherwise errors
 * Remarks		 : Executes at work-item level, meaning that only one work-item from work-group
                   executes this function at a time.
 */        
uint SNIP(TLLNode *pNodePool,	/* Hash-table + stack of free nodes */	
		  uint	   HeadIndex)	/* pointer to free-node in pNodePool */
{
    // private variables
    uint uiLLMNode = 0 ;/* 2:0->[f,r,m] and 31:3->next-node-index */
    uint uiLLNode = 0 ; /* next-node-index */
    uint pBits ; /* current-node-<f,r,m> bits, only 3 bits are using */;
    uint uiLLMNextNode = 0 ;/* 2:0->[f,r,m] and 31:3->next2next-node-index */
    uint uiLLNextNode = 0 ; /* next2next-node-index */
    uint nBits = 0 ; /* next-node-<f,r,m> bits, only 3 bits are using */;
    bool status = false ;

    uiLLMNode = pNodePool[HeadIndex].uiNext ; // next-node-index + <f,r,m>
    uiLLNode = GET_PTR(uiLLMNode) ; // next-node-index value, 31:3 bits
    pBits = GET_BITS(uiLLMNode) ; 
    if( ((!IS_SET(pBits,2)) & (!IS_SET(pBits,1))) && uiLLNode ) { // both retain,free has to be cleared
        uiLLMNextNode = pNodePool[uiLLNode].uiNext ; // 2:0->[f,r,m] and 31:3->next2next-node-index 
        uiLLNextNode = GET_PTR(uiLLMNextNode) ;
        nBits = GET_BITS(uiLLMNextNode) ; 
        if( (!IS_SET(nBits,2)) & (IS_SET(nBits,0)) ) { // <f,r,m> = <0,x,1>
            uiLLMNextNode = SET_PTR(uiLLNextNode) | pBits ; // make current-node next pointer as next2next node index
            atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[HeadIndex].uiNext));
            status = atomic_compare_exchange_strong(pChgPtr, &uiLLMNode, uiLLMNextNode) ;
            if(status == true) {
                uiFree(pNodePool, HeadIndex, uiLLNode) ;
                return HTS_SNIP_SUCCESS ; /* successfully deleted next marked node */
            }
            return HTS_SNIP_FAILED ; /* failed in deleting */						
        }
        return HTS_NEXT_UNMARKED ; /* next-node is not marked */
    } 
    return HTS_INVALID_START ; /* provided node is not valid one */
} // end-of-SNIP

/************************* CLEAN *****************************/
/* Function name : CLEAN
 * Arguments	 : HeadIndex, &next-node-index
 * Description	 : This function physically deletes node that are logically deleted.
                    each logically deleted node is snipped using SNIP function
 * Remarks		 : Runs at work-item-level
 */
uint CLEAN(TLLNode *pNodePool, /* Hash-table + stack of free nodes */
		   uint HeadIndex, /* pointer to free-node in pNodePool */
		   int *nextNodeIndex)
{
    uint uiLLNode = 0 ;
    uint status = 0, pBits = 0 ;
    uint uiLLMNextNode = 0 , uiLLNextNode = 0 ;

    uiLLNode = HeadIndex ;
    *nextNodeIndex = NULL ;

    while(true) {
        uiLLMNextNode = pNodePool[uiLLNode].uiNext ;
        uiLLNextNode  = GET_PTR(uiLLMNextNode) ;
        pBits = GET_BITS(uiLLMNextNode) ;
        status = SNIP(pNodePool, uiLLNode) ; // do it one-by-one
        if(status == HTS_NEXT_UNMARKED) {
            *nextNodeIndex = uiLLNextNode ; /* move to next node */
            return status ;
        } else if(status == HTS_INVALID_START) {
            if(( !IS_SET(pBits,2)) & (IS_SET(pBits,1))) {
                uiLLNode = uiLLNextNode ;
            } else { 
                if(uiLLNode == HeadIndex) {
                    return HTS_INVALID_START ;
                } else uiLLNode = HeadIndex ;
            } 
        } // end-of-first-if		
    } // end-of-while
} // end-of-CLEAN

/************************* REPLACE *****************************/ 
/* Function name : REPLACE
 * Arguments	 : pNodePool, HeadIndex, newIndex
 * Description	 : This function replaces HeadIndex with newIndex, if required 
                   replaces consequent nodes of newIndexa also.
 * Remarks		 : work-item-level function
 */
uint REPLACE(TLLNode *pNodePool, /* Hash-table + stack of free nodes */
			 uint HeadIndex, /* pointer to free-node in pNodePool */
			 uint uicurrentIndex,
			 int newIndex)
{
    // private variables
    uint uiLLMNode = 0 ;/* 2:0->[f,r,m] and 31:3->next-node-index */
    uint uiLLNode = 0 ; /* next-node-index */
    uint pBits ; /* current-node-<f,r,m> bits, only 3 bits are using */;
    uint uiLLMNextNode = 0 ;/* 2:0->[f,r,m] and 31:3->next2next-node-index */
    uint uiLLNextNode = 0 ; /* next2next-node-index */
    uint nBits = 0 ; /* next-node-<f,r,m> bits, only 3 bits are using */;
    bool status = false ;
    uint uiMtempNode = 0, uitempNode = 0 ;
    uint uitempIndex = 0 ;

    uiLLMNode = pNodePool[uicurrentIndex].uiNext ; // next-node-index + <f,r,m>
    uiLLNode = GET_PTR(uiLLMNode) ; // next-node-index value, 31:3 bits
    pBits = GET_BITS(uiLLMNode) ;
	
    // According to algorithm <f,r,m> has to be clear (i.e., 0)
    if( ((!IS_SET(pBits,2)) & (!IS_SET(pBits,1)) & (!IS_SET(pBits,0))) && uiLLNode ) { // all 3 bits should be clear

        // Get next-node details
        uiLLMNextNode = pNodePool[uiLLNode].uiNext ; // 2:0->[f,r,m] and 31:3->next2next-node-index 
        uiLLNextNode = GET_PTR(uiLLMNextNode) ;
        nBits = GET_BITS(uiLLMNextNode) ; 

        if( (!IS_SET(nBits,2)) & (!IS_SET(nBits,0)) ) { // <f,r,m> = <0,x,0>

            atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[uiLLNode].uiNext));
            uiMtempNode = SET_PTR(uiLLNextNode) | 0x2 ;  /* set retain bit */
            status = atomic_compare_exchange_strong(pChgPtr, &uiLLMNextNode, uiMtempNode) ;
            if(status == false) return HTS_RETAIN_FAILED ; /* failed */

        } else return HTS_INVALID_NBITS ;

        // who knows that new-node next pointer is 0 or not ??
        uitempIndex = newIndex ;
        uiMtempNode = pNodePool[uitempIndex].uiNext ;
        uitempNode = GET_PTR(uiMtempNode) ;
        while(uitempNode) {
            uitempIndex = uitempNode ;
            uiMtempNode = pNodePool[uitempIndex].uiNext ;
            uitempNode = GET_PTR(uiMtempNode) ;
        } ; // trace the path until end-node occurs
        
        pNodePool[uitempIndex].uiNext = SET_PTR(uiLLNextNode) | 0x0 ;
        atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[uicurrentIndex].uiNext));
        uiLLMNextNode = SET_PTR(newIndex) | pBits ;
        status = atomic_compare_exchange_strong(pChgPtr, &uiLLMNode, uiLLMNextNode) ;
        if(status == true) {
            uiFree(pNodePool, HeadIndex, uiLLNode) ; 
            return HTS_REPLACE_SUCCESS ;
        }
        else {
            atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[uiLLNode].uiNext));
            uiMtempNode = SET_PTR(uiLLNextNode) | 0x2 ;  /* set retain bit */
            status = atomic_compare_exchange_strong(pChgPtr, &uiLLMNextNode, uiMtempNode) ;
            return HTS_REPLACE_FAILED ;
        }
    } 
    /*else if( !uiLLNode ) { // boundary condition 
		// who knows that new-node next pointer is 0 or not ??
        uitempIndex = newIndex ;
        uiMtempNode = pNodePool[uitempIndex].uiNext ;
        uitempNode = GET_PTR(uiMtempNode) ;
        while(uitempNode) {
            uitempIndex = uitempNode ;
            uiMtempNode = pNodePool[uitempIndex].uiNext ;
            uitempNode = GET_PTR(uiMtempNode) ;
        } ; // trace the path until end-node occurs
        
        pNodePool[uitempIndex].uiNext = SET_PTR(uiLLNode) ;
        atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[uicurrentIndex].uiNext));
        uiLLMNextNode = SET_PTR(newIndex) | pBits ;
        status = atomic_compare_exchange_strong(pChgPtr, &uiLLMNode, uiLLMNextNode) ;
        if(status == true) return HTS_REPLACE_SUCCESS ;
        else return HTS_REPLACE_FAILED ;

	} */else return HTS_INVALID_PBITS ;

} // end-of-REPLACE

/************************* WINDOW *****************************/
 /* Function name : WINDOW
  * Arguments	  : key, prevNodeIndex, nextNodeIndex, &index
  * Description	  : This function returns index of node that a given key can be inserted
  * Remarks		  : work-group-level function
  */
uint WINDOW(TLLNode *pNodePool,	/* Hash-table + stack of free nodes */
			uint key,
			int prevNodeIndex, 
            int nextNodeIndex,
			uint *index)
{
    uint uiLLMNode = 0, uiLLNode = 0 ;
    uint uiLLMNextNode = 0, uiLLNextNode = 0 ;
    uint status = 0 ;
    uint pBits = 0, nBits = 0 ;
    uint pMaxVal = 0 ;

    uint lid = get_local_id(0) ;

    if(prevNodeIndex != (int)NULL) {
        uiLLMNode = pNodePool[prevNodeIndex].uiNext ;
        uiLLNode = GET_PTR(uiLLMNode) ;
        pBits = GET_BITS(uiLLMNode) ;
        if( IS_SET(pBits,2) || IS_SET(pBits,0) ) 
            return HTS_INVALID_PREV_BITS ;
        if( uiLLNode != nextNodeIndex) 
            return HTS_INVALID_NEXT_REF ;
        uint uiVal = pNodePool[prevNodeIndex].pE[lid];
        pMaxVal = work_group_reduce_max(uiVal) ;
        if(key <= pMaxVal) return HTS_WINDOW_NOT_FOUND ;							
    }

    uiLLMNextNode = pNodePool[nextNodeIndex].uiNext ;
    uiLLNextNode = GET_PTR(uiLLMNextNode) ;
    nBits = GET_BITS(uiLLMNextNode) ;
    if( IS_SET(nBits,0) || IS_SET(nBits,2) ) {
        return HTS_INVALID_NEXT_BITS ;
    }
    uint uiVal = pNodePool[nextNodeIndex].pE[lid];

    pMaxVal = work_group_reduce_max(uiVal) ;
    if(key <= pMaxVal) { 
        if(key == uiVal) {
            *index = lid + 1 ;
            status = HTS_KEY_FOUND ;
        } else {
            *index = 0 ;
             status = HTS_WINDOW_FOUND ;
        }
        *index = work_group_reduce_max(*index);
        if(*index) { // broadcast key 
            (*index) = work_group_broadcast(*index, (*index)-1) ;
            status = work_group_broadcast(status,(*index)-1) ;
        }
        work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
        return status ;
    } 
	return HTS_WINDOW_NOT_FOUND ;
} // end-of WINDOW

/* This function returns index of a key */
uint work_group_index(TLLNode *pNodePool,uint uiLLNode, uint key) {
     __local int keyIndex ; // cannot initialize local variable, getting error in compilation
     uint lid = get_local_id(0);
     if(lid == 0) keyIndex = 0 ; 
     work_group_barrier(CLK_LOCAL_MEM_FENCE);
     uint uiVal = pNodePool[uiLLNode].pE[lid] ;
     if(uiVal==key) { // thread conflicts does not occur because no duplication of keys
        keyIndex = (int)lid+1 ;
     }
     work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
     return keyIndex ;
}// end of work_group_index
/* This function returns minimum index of an EMPTY-KEY */
uint min_work_group_index(uint uiVal, uint maxIndex) {
	uint keyIndex = 0 ;
	uint lid = get_local_id(0) ;

	if(uiVal == EMPTY_KEY) {
		keyIndex = lid+1 ;
	} else keyIndex = maxIndex + 1;  // if 0, then below line always returns ZERO only
	keyIndex = work_group_reduce_min(keyIndex) ;
	keyIndex = work_group_broadcast(keyIndex, keyIndex-1) ;
	work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
	return keyIndex ;
}// end-of min_work_group_index

/************************* CLONE *****************************/
/* Function name : CLONE
 * Arguments	 : LLNode, key
 * Description	 : Clones a node and adds a key to it, if does not exist,
 *				   If a key exists, it deletes
 * Remarks		 : work-group level
 */  
 int CLONE(TLLNode *pNodePool,	/* Hash-table + stack of free nodes */
			uint uiHeadIndex,	/* pointer to free-node in pNodePool */
			uint nextNodeIndex,
			uint key)			/* key to inserted or deleted */
{
    int uiNewPtr = 0, uiprevNewPtr = 0 ;
    uint uiVal = 0, keyIndex = 0 ;
    uint lid = get_local_id(0) ;
    uint maxIndex = get_local_size(0) ; /* instead of 0, writing big value */

    if(lid == 0) {
        uiNewPtr = uiAlloc(pNodePool, uiHeadIndex) ;
    }
    work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
    uiNewPtr  = work_group_broadcast(uiNewPtr,0);
	// take care of boundary conditions here : if uiNewPtr is NULL
	if(uiNewPtr == NULL) {
		return uiNewPtr ;
	}
	pNodePool[uiNewPtr].uiNext = pNodePool[nextNodeIndex].uiNext ;

	keyIndex = work_group_index(pNodePool, nextNodeIndex, key) ;
	uiVal = pNodePool[nextNodeIndex].pE[lid] ;
 
    if(keyIndex == 0) { /* key-does-not-exist-insert-it */
        keyIndex = min_work_group_index(uiVal, maxIndex) ;
        // got minimum index of EMPTY_KEY 
        if(keyIndex <= maxIndex ) { // take care of extreme condition
            pNodePool[uiNewPtr].pE[lid] = pNodePool[nextNodeIndex].pE[lid];
            if(lid == (keyIndex-1)) {
                pNodePool[uiNewPtr].pE[lid] = key ;
            }
            work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
			return uiNewPtr ;	
        } else { /* node-is-full-create-new-node-for-storing-new-key*/
            if(lid == 0) {
                uiprevNewPtr = uiAlloc(pNodePool, uiHeadIndex) ;
            }
            work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);
            uiprevNewPtr = work_group_broadcast(uiprevNewPtr,0);
			// take care of boundary conditions here : if uiNewPtr is NULL
			if(uiNewPtr == NULL) {
				return uiNewPtr ;
			}
			// take care of boundary conditions here : if uiNewPtr is null		
            pNodePool[uiprevNewPtr].uiNext = SET_PTR(uiNewPtr) ;
            uiVal = pNodePool[nextNodeIndex].pE[lid] ;
            if(uiVal < key) /* values-lessthan-key-store-in-uiprevNewPtr */
                pNodePool[uiprevNewPtr].pE[lid] = uiVal ;
            else if(uiVal > key) /* values-greaterthan-key-store-in-uiNewPtr */
                pNodePool[uiNewPtr].pE[lid] = uiVal ;
            uiVal = pNodePool[uiNewPtr].pE[lid] ;
            // get-minimum-free-index */
            keyIndex = min_work_group_index(uiVal, maxIndex) ;
            if(lid == (keyIndex-1)) 
            pNodePool[uiNewPtr].pE[lid] = key ; /* store-actual-key */
            return uiprevNewPtr ;
          }
    } else { /* Key-exists-delete-it*/
        pNodePool[uiNewPtr].pE[lid] = pNodePool[nextNodeIndex].pE[lid];
        pNodePool[uiNewPtr].pE[(keyIndex-1)] = EMPTY_KEY ;
        return uiNewPtr ;
    }
    return NULL ;
 } // end of CLONE

/************************* TRAVERSE *****************************/
 /* Function name : TRAVERSE
  * Arguments	  : Previous-node-index, next-node-index
  * Descrition	  : traverses list from previous node to next node,
  *					executes at work-item level
  */
int TRAVERSE( TLLNode *pNodePool,	/* Hash-table + stack of free nodes */
			   int uiprevNodeIndex,
			   int uinextNodeIndex )
{
	uint uiPPtr = uiprevNodeIndex ;
	uint uiLLNode = 0, uiLLMNode = 0, pBits = 0 ;

	if(uiprevNodeIndex == (int)NULL) return uiprevNodeIndex ;

	uiLLMNode = pNodePool[uiPPtr].uiNext ;
	uiLLNode = GET_PTR(uiLLMNode) ; // get-node-inext
	pBits = GET_BITS(uiLLMNode) ; // get-bits <f,r,m>
	
	while( !IS_SET(pBits,2) ) {

		if( uiLLNode != uinextNodeIndex) {
			uiPPtr = uiLLNode ;
			// unpack next reference of uiPPtr
			uiLLMNode = pNodePool[uiPPtr].uiNext ;
			uiLLNode = GET_PTR(uiLLMNode) ;
			pBits = GET_BITS(uiLLMNode) ;
		} else {
			return uiPPtr ;
		} // end-of-if
	} // end-of-while
	return NULL ;
} // end of TRAVERSE

/*
** uiHashFunction:
** hash function to map key to hash table index.
*/
uint uiHashFunction(uint uiKey) 
{
    return uiKey & OCL_HASH_TABLE_MASK ;
}

/************************* bFind *****************************/
/* Function name	: bFind
 * Inputs			: NodePool, key
 * Outputs			: previousNode, NextNode, Index of key(if found)
 * Description		: work-group level function. searches for Key in the set, 
 *					  if key found return index of it, otherwise prevNode and nextNode
					  where the key can be inserted.	
 */
uint bFind( TLLNode* pNodePool,		/* Hash-table + stack of free nodes */		
		   uint     Key,			/* Key to be find */
           int*	prevNodeIndex,  /* if key not found, return window */
           int*	nextNodeIndex,
           uint*    Index)			/* if key found, return its index */
{
    int pRef = 0 ;
	int nRef = NULL ;
    uint status = 0, status2 = 0 ;
    bool bNodeFound = false ;
	uint uiLLMNextNode = 0, uiLLNextNode ;

    // get the thread id
    uint lid = get_local_id(0) ;

    // get the starting node
    uint uiPPtr  = uiHashFunction(Key) ;

    *prevNodeIndex = NULL ;
	*nextNodeIndex = NULL ;
    *Index = 0 ;

    pRef = NULL ;
    nRef = uiPPtr ;

	status = HTS_WINDOW_NOT_FOUND ;

    while(status == HTS_WINDOW_NOT_FOUND) {
		status = WINDOW(pNodePool, Key, pRef, nRef, Index) ;
		// window not found !!!
		if(status == HTS_WINDOW_NOT_FOUND) {
			// TODO : check on priority
			uiLLMNextNode = pNodePool[nRef].uiNext ;
			uiLLNextNode = GET_PTR(uiLLMNextNode) ;
			if(!uiLLNextNode) {// if next node doesnot exist, wht is d use of CLEAN
				status = HTS_WINDOW_FOUND ;
			}
			else { 
				pRef = nRef ;
				if(lid == 0) {
					status2 = CLEAN(pNodePool, pRef, &nRef) ;
					if(status2==HTS_INVALID_START)
						status = HTS_WINDOW_NOT_FOUND;
				}
			}
			work_group_barrier(CLK_GLOBAL_MEM_FENCE);
			status = work_group_broadcast(status,0) ;
			nRef = work_group_broadcast(nRef,0) ; 
		}
		// Handling boundary conditions --> TODO : check
		if( nRef == (int)NULL ) { // if this condition is not there, will endup in infinite loop
			status = HTS_WINDOW_FOUND ; 
		}
    } // end-of-while*/

    if(status == HTS_WINDOW_FOUND) { /* key not found, but window found */
        *prevNodeIndex = pRef ;
        *nextNodeIndex = nRef ;
        return status ;
    }
    if(status == HTS_KEY_FOUND) { /* key found */
		*prevNodeIndex = pRef ;
		*nextNodeIndex = nRef ;
        return status ;
    }

    return status ;
} // end of bFind()

/************************* bAdd *****************************/
/* Function name : bAdd
 * Arguments	 : NodePool, key, pointer-to-free-node-pool 
 * Description	 : work-group level function. searches for uiKey in the set, 
 *				   if key exists returns false, if does not exist inserts it 
 *				   in the required position.
 * returns true if uiKey is successfully inserted otherwise returns false
 */
bool bAdd( TLLNode*  pNodePool,		/* Hash-table + stack of free nodes */
		   uint		 uiKey,			/*Actual key to be inserted */
           uint      uiHeadIndex )	/* pointer to free-node in pNodePool */
{
	uint status  = 0, status2 = 0 ;
	int uiLLNode = NULL ;
	int prevNodeIndex = NULL, nextNodeIndex = NULL ;
	uint Index = 0 ;
	int cloneNodeIndex = 0 ;
	uint uiLLMNode = 0, uiLLNOde = 0, uiMNextNode = 0 ;
	uint pBits = 0 ;
	uint lid = get_local_id(0) ;

	// check whether key exists or not, if not exists return window_found
	// prevNodeIndex = nextNodEIndex in case of first element insertion
    status = bFind( pNodePool, uiKey, &prevNodeIndex, &nextNodeIndex, &Index ) ;

    if(status == HTS_WINDOW_FOUND) { // if key_does-not-found 
		
        if(lid == 0) {  // Go to next node 
            uiLLNode = TRAVERSE(pNodePool, prevNodeIndex, nextNodeIndex) ;
			//printf("uiLLNode%dprevNodeIndex%dnextNodeIndex%d", uiLLNode, prevNodeIndex, nextNodeIndex) ;
        }
		// wait for all work-items to reach this, point (especially WI0)
        work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE) ;
		// Let other work-items in a work-group know next-node-index !!!
        uiLLNode = work_group_broadcast(uiLLNode, 0) ;

        if( uiLLNode != (int)NULL) { // if next-node exists ?? 
			// clone the node by inserting key and return new node index
		    cloneNodeIndex = CLONE(pNodePool, uiHeadIndex, nextNodeIndex, uiKey) ;
			if(cloneNodeIndex == NULL) return false ; // Boundary conditions
			// Delete node which was there before inserting key and free it
            if(lid == 0) {
                status = REPLACE(pNodePool, uiHeadIndex, uiLLNode, cloneNodeIndex) ;
            }
			// wait for all work-items to reach this, point (especially WI0)
            work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE) ;
			// Let other work-items in a work-group know status value !!!
            status = work_group_broadcast(status, 0) ;

            if(status == HTS_REPLACE_SUCCESS) return true ; // SUCCESSFULLY INSERTED KEY 
        }//end-of if(uiLLNode != (int)NULL)
		else {
			// clone the node by inserting key and return new node index
		    cloneNodeIndex = CLONE(pNodePool, uiHeadIndex, nextNodeIndex, uiKey) ;
			if(cloneNodeIndex == NULL) return false ; // Boundary conditions
			uiLLMNode = pNodePool[nextNodeIndex].uiNext ;
			uiLLNode = GET_PTR(uiLLMNode) ;
			pBits = GET_BITS(uiLLMNode) ;
			if(lid == 0) {
				pNodePool[cloneNodeIndex].uiNext = pNodePool[nextNodeIndex].uiNext ;
				atomic_uint* pChgPtr =
                (atomic_uint *)(&(pNodePool[nextNodeIndex].uiNext));
				uiMNextNode = SET_PTR(cloneNodeIndex) | pBits ;
				status2 = atomic_compare_exchange_strong(pChgPtr, &uiLLMNode, uiMNextNode) ;
				if(status2 == false) status = HTS_REPLACE_FAILED ;
				else status = HTS_REPLACE_SUCCESS ;
			}
			// wait for all work-items to reach this, point (especially WI0)
            work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE) ;
			// Let other work-items in a work-group know status value !!!
            status = work_group_broadcast(status, 0) ;

            if(status == HTS_REPLACE_SUCCESS) return true ; // SUCCESSFULLY INSERTED KEY 
		}

    }// end-of if(status == HTS_WINDOW_FOUND)
	return false ; /* FAILED IN INSERTING KEY */

} // end of bAdd()


/************************* bRemove *****************************/
/* Function name : bRemove
 * Arguments	 : NodePool, key, pointer-to-free-node-pool 
 * Description	 : work-group level function. searches for uiKey in the set, 
 *				   if key does not exist, returns false otherwise deletes it.
 * returns true if uiKey is successfully deleted otherwise returns false
 */
bool bRemove( TLLNode*  pNodePool,    /* Hash-table + stack of free nodes */
			  uint	    uiKey,	      /* Actual key to be deleted */
              uint      uiHeadIndex ) /* pointer to free-node in pNodePool */
{
	uint status  = 0 ;
	int uiLLNode = NULL ;
	uint prevNodeIndex = 0, nextNodeIndex = 0, Index = 0 ;
	int cloneNodeIndex = 0 ;
	uint lid = get_local_id(0) ;

	// check whether key exists or not, if exists return key_found
    status = bFind( pNodePool, uiKey, &prevNodeIndex, &nextNodeIndex, &Index ) ;

    if(status == HTS_KEY_FOUND) { /* if key_found */

        if(lid == 0) {  /* Go to next node */
            uiLLNode = TRAVERSE(pNodePool, prevNodeIndex, nextNodeIndex) ;
			//printf("uiLLNode%dprevNodeIndex%dnextNodeIndex%d", uiLLNode, prevNodeIndex, nextNodeIndex) ;
        }
		// wait for all work-items to reach this, point (especially WI0)
        work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE) ;
		// Let other work-items in a work-group know next-node-index !!!
        uiLLNode = work_group_broadcast(uiLLNode, 0) ;

        if( uiLLNode != (int)NULL) { /* if next-node exists ?? */

			// clone the node by deleting key and return new node
		    cloneNodeIndex = CLONE(pNodePool,uiHeadIndex,nextNodeIndex,uiKey) ;
			if(cloneNodeIndex == NULL) return false ; // Boundary conditions
			
			// Delete node which was there before deleting key and free it
            if(lid == 0) {
                status = REPLACE(pNodePool,uiHeadIndex, uiLLNode, cloneNodeIndex) ;
            }

			// wait for all work-items to reach this, point (especially WI0)
            work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE) ;
			// Let other work-items in a work-group know status value !!!
            status = work_group_broadcast(status, 0) ;

            if(status == HTS_REPLACE_SUCCESS) return true ; /* SUCCESSFULLY DELETED KEY */

        }//end-of if(uiLLNode != (int)NULL)

    }// end-of if(status == HTS_WINDOW_FOUND)

	return false ; /* FAILED IN DELETING KEY */
} // end of bRemove()


// kernel calling from host
__kernel void HTSTopKernel(__global void* pvOclReqQueue,
                           __global void* pvNodePool,
                           __global void* pvMiscData,
                           uint  uiReqCount)
    {
    uint uiFlags = 0 ;
    uint uiKey = 0 ;
    uint uiStatus = 0 ;
    uint uiType = 0 ;
    uint uiretVal = 0 ; 
    bool bReqStatus = false ;
    
    //get the svm data structures
    TQueuedRequest* pOclReqQueue = (TQueuedRequest *)pvOclReqQueue ;
    TLLNode*        pNodePool    = (TLLNode *)pvNodePool ; // hash table + free node pool
    TMiscData*      pMiscData    = (TMiscData*)pvMiscData ;
    uint            uiHeadIndex =  pMiscData->uiHeadIndex; // points to head of free nodes pool

    uint grid = get_group_id(0);
    uint lid  = get_local_id(0);

    if (grid < uiReqCount)
        {
        uiKey   = pOclReqQueue[grid].uiKey; 
        uiType  = pOclReqQueue[grid].uiType;
        uiFlags = pOclReqQueue[grid].uiFlags;

        uint uiNode = 0, uiIndex = 0 ;
        int prevNodeIndex = NULL, nextNodeIndex = NULL ;

        if(uiType == HTS_REQ_TYPE_FIND) 
            {
				/* when key_found, no need to worry about prevNodeIndex and nextNodeIndex */
				prevNodeIndex = 0; nextNodeIndex = 0 ; uiIndex = 0 ;
				uiretVal = bFind(pNodePool,
							uiKey,
							&prevNodeIndex,
							&nextNodeIndex, /* key location */
							&uiIndex);
                if(uiIndex) bReqStatus = true ;
				if(lid == 0) printf("uiIndex%d-prevNodeIndex%d-nextNodeIndex%d",uiIndex,prevNodeIndex, nextNodeIndex) ;
            }
        else if(uiType == HTS_REQ_TYPE_ADD)
            {
                bReqStatus = bAdd(pNodePool, 
								  uiKey, 
								  uiHeadIndex) ;
            }
        else if(uiType == HTS_REQ_TYPE_REMOVE)
            {
                bReqStatus = bRemove(pNodePool, 
								  uiKey, 
								  uiHeadIndex) ;
            }
        work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

        if(bReqStatus == true)
            uiIndex = 1;
        else
            uiIndex = 0;

        if(lid == 0)
            {
            uiFlags = SET_FLAG(uiFlags,HTS_REQ_COMPLETED);
            pOclReqQueue[grid].uiFlags  = uiFlags;
            pOclReqQueue[grid].uiStatus = uiIndex;
            }
        work_group_barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);

        } 
    }



